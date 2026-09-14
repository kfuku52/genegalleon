"""Anonymous public plant-trait adapters; source rows remain inspectable on disk."""
import hashlib
import io
import json
import os
import re
import tempfile
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import quote
from urllib.request import Request

import pandas as pd
from format_species_network import guarded_urlopen

DATABASES = ('bien', 'brot', 'cpt', 'algaetraits')
CPT_FILES = {
    'species': 38095881, 'taxonomy': 38217405, 'chemical': 38095812, 'hydraulic': 38095830,
    'morphometric': 38095836, 'photosynthetic': 38095848,
    'functional': 38095854, 'pathway': 42349308,
}
COLUMNS = ['species', 'trait_key', 'value', 'unit', 'source_id', 'source_reference', 'context']


def atomic_write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(data)
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


class Acquisition:
    def __init__(self, directory, timeout, opener=guarded_urlopen):
        self.directory = Path(directory)
        self.timeout = timeout
        self.opener = opener
        self.records = []
        self.decisions = []
        self.last_request = None

    def get(self, url, *, bien_empty=False):
        # The published BIEN backend limits an IP to 100 requests / 15 minutes.
        interval = 9.1 if bien_empty else 0.4
        if self.last_request is not None:
            time.sleep(max(0, interval - (time.monotonic() - self.last_request)))
        self.last_request = time.monotonic()
        try:
            with self.opener(Request(url, headers={'User-Agent': 'genegalleon-public-traits'}), timeout=self.timeout) as response:
                data = response.read()
        except HTTPError as exc:
            # BIEN's documented no-observations response is distinct from a missing route.
            if not bien_empty or exc.code != 404:
                raise
            data = exc.read()
            try:
                payload = json.loads(data)
            except ValueError:
                raise ValueError('BIEN endpoint unavailable (non-JSON 404)') from exc
            if not isinstance(payload, dict) or not str(payload.get('message', '')).startswith(('No data found for species ', 'No traits available for species ')):
                raise ValueError('Unexpected BIEN 404 response') from exc
            self.save(url, data)
            return b''
        self.save(url, data)
        return data

    def save(self, url, data):
        digest = hashlib.sha256(data).hexdigest()
        path = self.directory / 'raw' / (digest + '.bin')
        atomic_write(path, data)
        self.records.append({'url': url, 'sha256': digest, 'path': str(path),
                             'retrieved_at': datetime.now(timezone.utc).isoformat()})

    def json(self, url, **kwargs):
        data = self.get(url, **kwargs)
        return json.loads(data) if data.strip() else []

    def csv(self, url, encoding=None):
        data = self.get(url)
        # Chinese release CSVs use a legacy encoding; Latin-1 preserves bytes,
        # including the hydraulic header, without replacement characters.
        try:
            text = data.decode(encoding or 'utf-8-sig')
        except UnicodeDecodeError:
            text = data.decode('latin1')
        return pd.read_csv(io.StringIO(text), dtype=str, keep_default_na=False)


def exact_species(name):
    """Do not promote qualified names, morphospecies or higher taxa to species."""
    name = str(name).strip().replace('_', ' ')
    return name if re.fullmatch(r'[A-Z][a-z-]+ [a-z][a-z-]+', name) and name.split()[1] not in {'sp', 'spp', 'cf', 'aff'} else None


def context(row):
    return json.dumps(dict(row), ensure_ascii=False, sort_keys=True)


def brot(acquisition, config, targets):
    table = acquisition.csv(config.get('uri') or 'https://ndownloader.figshare.com/files/11194784')
    required = {'ID', 'Taxon', 'Trait', 'Data', 'Units', 'DataType', 'SourceID'}
    if not required.issubset(table.columns):
        raise ValueError('BROT source is missing required columns')
    rows = []
    for record in table.to_dict('records'):
        name = exact_species(record['Taxon'])
        if name not in targets:
            continue
        # BROT mixes quantitative and categorical representations of the same
        # trait. Never pool percentages, yes/no, or different physical units.
        key = '|'.join([record['Trait'], record['DataType'], record['Units']])
        rows.append([targets[name], key, record['Data'], record['Units'], record['ID'], record['SourceID'], context(record)])
    return rows


def cpt(acquisition, config, targets):
    base = config.get('uri') or 'https://ndownloader.figshare.com/files/'
    tables = {name: acquisition.csv(base.rstrip('/') + '/' + str(file_id), **({'encoding': 'gb18030'} if name == 'hydraulic' else {})) for name, file_id in CPT_FILES.items()}
    mapping = tables['species']
    required = {'SAMPLE ID', 'SPECIES ID', 'ACCEPTED GENUS', 'ACCEPTED SPECIES'}
    if not required.issubset(mapping.columns):
        raise ValueError('CPT species translations missing required columns')
    mapping = mapping.copy()
    mapping['species'] = (mapping['ACCEPTED GENUS'] + ' ' + mapping['ACCEPTED SPECIES']).map(exact_species)
    if mapping['SAMPLE ID'].duplicated().any() or mapping['SAMPLE ID'].eq('').any():
        raise ValueError('CPT sample mapping must have unique nonempty sample IDs')
    taxonomy = tables['taxonomy'].copy()
    if not {'SPECIES ID', 'ACCEPTED GENUS', 'ACCEPTED SPECIES'}.issubset(taxonomy.columns):
        raise ValueError('CPT taxonomy missing required columns')
    taxonomy['species'] = (taxonomy['ACCEPTED GENUS'] + ' ' + taxonomy['ACCEPTED SPECIES']).map(exact_species)
    rows = []
    for table_name, table in tables.items():
        if table_name in {'species', 'taxonomy'}:
            continue
        key = 'SPECIES ID' if table_name == 'pathway' else 'SAMPLE ID'
        if key not in table:
            raise ValueError('CPT table missing ' + key)
        lookup = (taxonomy if table_name == 'pathway' else mapping)[[key, 'species']].drop_duplicates()
        if lookup[key].duplicated().any():
            raise ValueError('CPT ambiguous species mapping')
        merged = table.merge(lookup, on=key, how='left', validate='many_to_one', indicator=True)
        if (merged['_merge'] != 'both').any():
            raise ValueError('CPT trait rows contain unmapped IDs')
        for record in merged.to_dict('records'):
            name = record['species']
            if name not in targets:
                continue
            # Published outlier flags apply to the original record. Keep raw
            # records, but exclude flagged observations from analytical values.
            flagged = str(record.get('flagged', '')).strip()
            if flagged not in {'', '0', 'NA'}:
                continue
            for column in table.columns:
                if column in {key, 'flagged'}:
                    continue
                value = record[column]
                rows.append([targets[name], table_name + ':' + column.strip(), value,
                             'see_CPT_v2_dictionary', table_name + ':' + record[key],
                             'https://doi.org/10.6084/m9.figshare.19448219', context(record)])
    return rows


def bien(acquisition, config, targets):
    base = config.get('uri') or 'https://mint-pheasant.nceas.ucsb.edu:5775/api/download/traits'
    response_format = config.get('response_format') or 'csv'
    if response_format not in {'csv', 'json'}:
        raise ValueError('BIEN response_format must be csv or json')
    rows = []
    for name, label in targets.items():
        url = base + '?species=' + quote(name, safe='')
        data = acquisition.get(url, bien_empty=True)
        if not data.strip():
            records = []
        elif response_format == 'json':
            payload = json.loads(data)
            if not isinstance(payload, dict) or not isinstance(payload.get('data'), list):
                raise ValueError('BIEN response must contain a data array')
            records = payload['data']
        else:
            table = pd.read_csv(io.StringIO(data.decode('utf-8-sig')), dtype=str, keep_default_na=False)
            if not {'trait_id', 'scrubbed_species_binomial', 'trait_name', 'trait_value', 'unit', 'access'}.issubset(table.columns):
                raise ValueError('BIEN CSV missing required columns')
            records = table.rename(columns={'trait_id': 'id'}).to_dict('records')
        for record in records:
            required = {'scrubbed_species_binomial', 'trait_name', 'trait_value', 'unit', 'id', 'access'}
            if not isinstance(record, dict) or not required.issubset(record):
                raise ValueError('BIEN record missing required fields')
            if record['scrubbed_species_binomial'] != name:
                raise ValueError('BIEN returned an unexpected species')
            if str(record['access']).strip().lower() != 'public':
                acquisition.decisions.append({'species': name, 'source_id': str(record['id']), 'status': 'excluded_nonpublic'})
                continue
            # Unit is part of the key to prevent scientifically invalid pooling.
            unit = str(record['unit'] or '')
            rows.append([label, record['trait_name'] + '|' + unit, record['trait_value'], unit,
                         str(record['id']), record.get('source_citation') or record.get('url_source', ''), context(record)])
    return rows


def algaetraits(acquisition, config, targets):
    base = (config.get('uri') or 'https://www.marinespecies.org/rest').rstrip('/')
    rows = []
    for name, label in targets.items():
        hits = acquisition.json(base + '/AphiaRecordsByName/' + quote(name, safe='') + '?like=false&marine_only=false')
        if not isinstance(hits, list):
            raise ValueError('WoRMS name lookup must return an array')
        matches = [hit for hit in hits if hit.get('scientificname') == name and hit.get('rank') == 'Species' and hit.get('status') == 'accepted']
        if len(matches) != 1:
            if hasattr(acquisition, 'decisions'):
                acquisition.decisions.append({'species': name, 'status': 'unresolved', 'exact_accepted_hits': len(matches)})
            continue
        aphia = matches[0]['AphiaID']
        attributes = acquisition.json(base + '/AphiaAttributesByAphiaID/' + str(aphia))
        if not isinstance(attributes, list):
            raise ValueError('WoRMS attributes must return an array')
        for record in attributes:
            if not {'AphiaID', 'measurementTypeID', 'measurementType', 'measurementValue', 'AphiaID_Inherited'}.issubset(record):
                raise ValueError('WoRMS attribute missing required fields')
            # Never inherit genus/other-species observations. Child attributes
            # describe the parent measurement (unit, locality, life stage).
            if str(record['AphiaID']) != str(aphia) or str(record['AphiaID_Inherited']) != str(aphia):
                continue
            qualifiers = record.get('children', [])
            key = str(record['measurementTypeID']) + '|' + record['measurementType']
            if qualifiers:
                # Preserve complex measurements as JSON, not unqualified numbers.
                value = context({'value': record['measurementValue'], 'qualifiers': qualifiers})
            else:
                value = record['measurementValue']
            rows.append([label, key, value, 'see_context', str(record.get('source_id', '')),
                         record.get('reference', ''), context(record)])
    return rows


def load_public_traits(database, config, species, downloads_dir, timeout, dry_run, opener=guarded_urlopen):
    if database not in DATABASES:
        raise ValueError('Unknown public plant database: ' + database)
    targets = {}
    for label in species:
        name = exact_species(label)
        if name:
            if name in targets and targets[name] != label:
                raise ValueError('Duplicate species labels for ' + name)
            targets[name] = label
    if dry_run:
        print('[dry-run] anonymous {} acquisition for {} exact species'.format(database, len(targets)))
        return None
    directory = Path(downloads_dir) / database
    acquisition = Acquisition(directory, timeout, opener)
    report = {'database': database, 'targets': list(species), 'status': 'running'}
    report['excluded_target_labels'] = [label for label in species if not exact_species(label)]
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f') + '-' + uuid.uuid4().hex
    try:
        rows = globals()[database](acquisition, config, targets)
        table = pd.DataFrame(rows, columns=COLUMNS).drop_duplicates()
        observed = table['value'].map(lambda value: value is not None and str(value).strip().lower() not in {'', 'na', 'nan', 'n/a', 'null', 'none', 'unknown'})
        report['missing_value_rows'] = int((~observed).sum())
        table = table.loc[observed].reset_index(drop=True)
        data = table.to_csv(sep='\t', index=False).encode()
        normalized = directory / 'runs' / (run_id + '.tsv')
        atomic_write(normalized, data)
        report.update(status='complete', rows=len(table), normalized=str(normalized),
                      normalized_sha256=hashlib.sha256(data).hexdigest())
        return table
    except Exception as exc:
        report.update(status='failed', error=str(exc))
        raise
    finally:
        report['requests'] = acquisition.records
        report['species_decisions'] = acquisition.decisions
        atomic_write(directory / 'runs' / (run_id + '.json'), json.dumps(report, indent=2).encode())

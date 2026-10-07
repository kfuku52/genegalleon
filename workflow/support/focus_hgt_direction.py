"""Filter modeled species branches using existing host-species taxonomy."""

import csv
import hashlib
from collections import Counter
from pathlib import Path

DIRECTION_CHOICES = ('any', 'non_arthropoda_to_insecta')
EVENT_FIELDS = ['direction_filter_status', 'direction_filter_reason',
                'direction_donor_arthropoda_status', 'direction_recipient_insecta_status']
MISSING = {'', '.', 'na', 'nan', 'none', 'null', 'unknown', 'unavailable'}


def key(value):
    return str(value).strip().replace(' ', '_')


def known(value):
    return str(value).strip().lower() not in MISSING


def read_taxonomy(path):
    raw = Path(path).read_bytes()
    import io
    reader = csv.DictReader(io.StringIO(raw.decode('utf-8-sig')), delimiter='\t')
    fields = reader.fieldnames or []
    required = {'species', 'resolution_status', 'domain', 'phylum', 'class'}
    if len(fields) != len(set(fields)) or not required <= set(fields):
        raise ValueError('Missing or duplicate species-taxonomy columns')
    aliases, species = {}, set()
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            raise ValueError('Malformed species-taxonomy row')
        name = key(row['species'])
        if not name or name in species:
            raise ValueError('Duplicate or empty taxonomy species identity')
        species.add(name)
        if row['resolution_status'] == 'resolved' and row['class'] == 'Insecta' and row['phylum'] != 'Arthropoda':
            raise ValueError('Insecta classification requires Arthropoda phylum: ' + name)
        names = {name}
        if row.get('tree_status') == 'mapped' and known(row.get('tree_tip', '')):
            names.add(key(row['tree_tip']))
        for alias in names:
            if alias in aliases:
                raise ValueError('Ambiguous taxonomy species/tree-tip alias: ' + alias)
            aliases[alias] = row
    return aliases, hashlib.sha256(raw).hexdigest()


def tip_state(row, group):
    if not row or row['resolution_status'] != 'resolved':
        return 'unknown'
    phylum, domain, klass = (row[name].strip() for name in ('phylum', 'domain', 'class'))
    if group == 'arthropoda':
        return 'within' if phylum == 'Arthropoda' else 'outside' if known(phylum) or domain in {'Bacteria', 'Archaea'} else 'unknown'
    if klass == 'Insecta':
        return 'within'
    return 'outside' if known(klass) or domain in {'Bacteria', 'Archaea'} or known(phylum) and phylum != 'Arthropoda' else 'unknown'


def classify_branches(nodes, taxonomy):
    rows, states = [], {}
    for name, tips in sorted(nodes.items()):
        result = dict(species_branch=name, clade_tip_labels='; '.join(tips), descendant_tip_count=len(tips))
        for group in ('arthropoda', 'insecta'):
            counts = Counter(tip_state(taxonomy.get(tip), group) for tip in tips)
            status = ('mixed' if counts['within'] and counts['outside'] else 'unknown' if counts['unknown'] or not tips
                      else 'within' if counts['within'] else 'outside')
            result[group + '_status'] = status
            for state in ('within', 'outside', 'unknown'):
                result[f'{group}_{state}_tip_count'] = counts[state]
            result[group + '_unknown_tip_labels'] = '; '.join(tip for tip in tips if tip_state(taxonomy.get(tip), group) == 'unknown')
        states[name] = result
        rows.append(result)
    return states, rows


def filter_events(events, nodes, taxonomy_path):
    """Require every donor tip outside Arthropoda and every recipient tip in Insecta."""
    taxonomy, source_sha = read_taxonomy(taxonomy_path)
    states, branches = classify_branches(nodes, taxonomy)
    selected, audit = [], []
    for event in events:
        donor, recipient = (key(event[f'generax_{side}_node']) for side in ('donor', 'recipient'))
        if event['generax_transfer'] != f"Y@{event['generax_donor_node']}@{event['generax_recipient_node']}":
            raise ValueError('Event transfer token disagrees with donor/recipient branch IDs')
        ds = states.get(donor, {}).get('arthropoda_status', 'unmapped')
        rs = states.get(recipient, {}).get('insecta_status', 'unmapped')
        reason = ('event_mapping_unresolved' if event.get('mapping_status', 'matched') != 'matched'
                  else 'donor_' + ds + '_arthropoda' if ds != 'outside'
                  else 'recipient_' + rs + '_insecta' if rs != 'within' else '')
        status = ('passed' if not reason else 'withheld' if event.get('mapping_status', 'matched') != 'matched'
                  or ds in {'unknown', 'mixed', 'unmapped'} or rs in {'unknown', 'mixed', 'unmapped'} else 'excluded_direction')
        row = dict(event, direction_filter_status=status,
                   direction_filter_reason=reason or 'all_donor_tips_outside_arthropoda_all_recipient_tips_in_insecta',
                   direction_donor_arthropoda_status=ds, direction_recipient_insecta_status=rs)
        audit.append(row)
        if status == 'passed':
            selected.append(row)
    return selected, audit, branches, source_sha

"""Independent transcript evidence from other accessions of a target species."""
import contextlib
import csv
import io
import re
import zipfile
from collections import defaultdict
from itertools import pairwise
from pathlib import Path

from cdskit.sequence_selection import fasta_records

try:
    from fasta_sequence_store import open_text
    from refinement_native_splits import collinear_gap
    from rescue_coding_paths import OwnershipIndex
except ImportError:
    from .fasta_sequence_store import open_text
    from .refinement_native_splits import collinear_gap
    from .rescue_coding_paths import OwnershipIndex

ROLES = ('gtf', 'protein')
FIELDS = ('species', 'target_accession', 'evidence_accession', *ROLES, 'rna_sources')
TRANSCRIPT = re.compile(r'(?:^|;\s*)transcript_id\s+"([^"]+)"')
GENE = re.compile(r'(?:^|;\s*)gene_id\s+"([^"]+)"')


def read_manifest(path, species):
    """Bind same-species accession evidence to target IDs, with no locus list."""
    entries = defaultdict(list)
    paths = set()
    accessions = set()
    targets = {}
    with path.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames or set(reader.fieldnames) != set(FIELDS) or len(reader.fieldnames) != len(FIELDS):
            raise ValueError('Accession RNA manifest requires columns: ' + ', '.join(FIELDS))
        for raw in reader:
            if None in raw or any(not value or not value.strip() for value in raw.values()):
                raise ValueError(f'Incomplete accession RNA row at line {reader.line_num}')
            row = {k: v.strip() for k, v in raw.items()}
            name, target, accession = (row[k] for k in FIELDS[:3])
            key = name, accession.casefold()
            if (name not in species or target.casefold() == key[1] or key in accessions
                    or targets.setdefault(name, target.casefold()) != target.casefold()):
                raise ValueError(f'Invalid or repeated accession RNA source at line {reader.line_num}')
            accessions.add(key)
            row['rna_sources'] = sorted({label.strip() for label in row['rna_sources'].split(',') if label.strip()})
            if not row['rna_sources']:
                raise ValueError(f'No evidence labels at line {reader.line_num}')
            for role in ROLES:
                resource = (path.parent / row[role]).resolve(strict=True)
                if not resource.is_file() or resource in paths:
                    raise ValueError('Invalid or reused accession RNA input: ' + str(resource))
                row[role] = str(resource)
                paths.add(resource)
            entries[name].append({k: v for k, v in row.items() if k != 'species'})
    if not entries:
        raise ValueError('Empty accession RNA manifest')
    return {name: sorted(rows, key=lambda row: row['evidence_accession']) for name, rows in entries.items()}


@contextlib.contextmanager
def open_resource(path, *, newline=None):
    """Stream GTF/FASTA files or a ZIP containing exactly one such text file."""
    if path.suffix.lower() != '.zip':
        with open_text(path, newline=newline) as handle:
            yield handle
        return
    with zipfile.ZipFile(path) as archive:
        entries = [entry for entry in archive.infolist() if not entry.is_dir()
                   and not entry.filename.startswith('__MACOSX/')
                   and Path(entry.filename).suffix.lower() in {'.gtf', '.fa', '.fasta', '.faa', '.pep'}]
        if len(entries) != 1:
            raise ValueError('Expected exactly one GTF or FASTA ZIP member: ' + str(path))
        with io.TextIOWrapper(archive.open(entries[0]), encoding='utf-8', newline=newline) as handle:
            yield handle


def exons(path, *, contigs=(), transcript_ids=()):
    contig_prefixes = tuple(seqid + '\t' for seqid in contigs)
    # A small set of exact quoted IDs avoids splitting and parsing most lines.
    # The attribute parser below still decides which ID owns each exon.
    match = (re.compile(r'"(?:' + '|'.join(re.escape(i) for i in sorted(transcript_ids)) + r')"').search
             if 0 < len(transcript_ids) <= 64 else None)
    with open_resource(path) as handle:
        for line in handle:
            if (contig_prefixes and not line.startswith(contig_prefixes)) or (match and not match(line)):
                continue
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t', 8)
            if len(fields) != 9:
                raise ValueError('Invalid accession RNA GTF row: ' + str(path))
            if fields[2] == 'exon':
                yield fields


def matching_transcripts(path, proteins):
    hits = defaultdict(set)
    identifiers = set()
    for identifier, _, protein in fasta_records(path, opener=open_resource):
        protein = protein.upper().removesuffix('*')
        if identifier in identifiers:
            raise ValueError('Duplicate accession RNA protein ID: ' + identifier)
        identifiers.add(identifier)
        if protein in proteins:
            hits[protein].add(identifier)
    return hits


def transcript_paths(path, identifiers):
    models = {}
    for seqid, source, _, start, end, _, strand, _, attrs in exons(path, transcript_ids=identifiers):
        transcript = TRANSCRIPT.search(attrs)
        if transcript is None or transcript[1] not in identifiers:
            continue
        gene = GENE.search(attrs)
        start, end = int(start) - 1, int(end)
        if gene is None or strand not in {'+', '-'} or not 0 <= start < end:
            raise ValueError('Invalid accession RNA exon geometry: ' + transcript[1])
        t = models.setdefault(transcript[1], {'locus': (gene[1], seqid, strand), 'blocks': set(), 'labels': set()})
        if t['locus'] != (gene[1], seqid, strand):
            raise ValueError('Conflicting accession RNA transcript ID: ' + transcript[1])
        t['blocks'].add((start, end))
        t['labels'].add(source)
    for model in models.values():
        blocks = model['blocks'] = sorted(model['blocks'])
        if any(end > start for (_, end), (start, _) in pairwise(blocks)):
            raise ValueError('Overlapping accession RNA exons: ' + model['locus'][0])
        model['position'] = (*model['locus'], blocks[0][0], blocks[-1][1])
        model['exonic_bases'] = sum(end - start for start, end in blocks)
    return models


def single_rna_locus(ids, models, labels, protein):
    if not ids or any(identifier not in models for identifier in ids):
        return []
    if len({models[identifier]['locus'] for identifier in ids}) != 1:
        return []
    return [(identifier, model) for identifier in sorted(ids) if (model := models[identifier])['labels'] <= labels
            and model['exonic_bases'] >= 3 * len(protein)]


def split_pair(left, right):
    """Require two distinct same-strand transcripts in transcriptional order."""
    for a_id, a in left:
        for b_id, b in right:
            gap = collinear_gap(a['position'], b['position'])
            if gap is not None and gap <= 10_000:
                return (a_id, a), (b_id, b)
    return None


def spanning_transcripts(path, pairs, labels):
    """Find RNA-backed exons of the same transcript on both sides."""
    features = []
    for key, pair in pairs.items():
        for side, (_, model) in enumerate(pair):
            _, seqid, strand, start, end = model['position']
            blocks = model['blocks']
            features.append({'gene_id': (key, side), 'seqid': seqid, 'strand': strand,
                             'start': start, 'end': end, 'cds': blocks})
    index = OwnershipIndex(features)
    bounds = {seqid: (starts[0], maxima[-1]) for seqid, (starts, _, maxima) in index.tracks.items()}
    seen = {}
    blocked = set()
    for seqid, source, _, start, end, _, strand, _, attrs in exons(path, contigs=bounds):
        if source not in labels:
            continue
        a, b = int(start) - 1, int(end)
        left, right = bounds[seqid]
        if b <= left or a >= right:
            continue
        if not 0 <= a < b:
            raise ValueError('Invalid accession RNA exon geometry at ' + seqid)
        overlap = index.overlapping({'seqid': seqid, 'strand': strand, 'start': a, 'end': b, 'cds': ((a, b),)})
        if not overlap:
            continue
        transcript, gene = TRANSCRIPT.search(attrs), GENE.search(attrs)
        if transcript is None or gene is None:
            raise ValueError('Accession RNA GTF exon lacks transcript_id or gene_id')
        for row in overlap:
            key, side = row['gene_id']
            if key in blocked:
                continue
            identity = key, seqid, strand, gene[1], transcript[1]
            seen[identity] = seen.get(identity, 0) | (1 << side)
            if seen[identity] == 3:
                blocked.add(key)
    return blocked


def accession_support(entries, proposals, eligible):
    """Match validated ORFs exactly, without conflating accessions with donors."""
    if not entries or not eligible:
        return {}
    owners = defaultdict(set)
    for candidate in proposals.values():
        owners[candidate['proteins']].add(candidate['owner'])
    pending = {key: proposals[key] for key in eligible if len(owners[proposals[key]['proteins']]) == 1}
    fused_proteins = {key: c['source_protein'].upper().removesuffix('*') for key, c in pending.items()}
    proteins = set(fused_proteins.values()) | {protein for c in pending.values() for protein in c['proteins']}
    if not proteins:
        return {}
    supported, opposed = defaultdict(list), set()
    for entry in entries:
        gtf = Path(entry['gtf'])
        labels = set(entry['rna_sources'])
        hits = matching_transcripts(Path(entry['protein']), proteins)
        models = transcript_paths(gtf, set().union(*hits.values())) if hits else {}
        pairs = {}
        for key, candidate in pending.items():
            fused = fused_proteins[key]
            if any(identifier in models and models[identifier]['labels'] <= labels
                   and models[identifier]['exonic_bases'] >= 3 * len(fused) for identifier in hits[fused]):
                opposed.add(key)
                continue
            a, b = (single_rna_locus(hits[p], models, labels, p) for p in candidate['proteins'])
            if pair := split_pair(a, b):
                pairs[key] = pair
        blocked = spanning_transcripts(gtf, pairs, labels) if pairs else set()
        opposed.update(blocked)
        for key, pair in pairs.items():
            if key not in blocked:
                supported[key].append({
                    'target_accession': entry['target_accession'], 'evidence_accession': entry['evidence_accession'],
                    'transcripts': [identifier for identifier, _ in pair],
                    'genes': [model['locus'][0] for _, model in pair],
                    'rna_sources': [sorted(model['labels']) for _, model in pair],
                })
    return {key: matches for key, matches in supported.items() if key not in opposed}

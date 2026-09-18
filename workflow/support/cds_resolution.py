#!/usr/bin/env python3
"""Resolve CDS candidates without changing source FASTA/GFF files.

The output is an analysis input, not a corrected biological annotation. Every
selection (including exclusion) retains source and translation evidence.
"""
from __future__ import annotations

import argparse
import fcntl
import gzip
import hashlib
import importlib.metadata
import json
import os
import re
import shutil
import tempfile
from pathlib import Path
from urllib.parse import unquote

import pandas as pd
from Bio.Data import CodonTable
from Bio.Seq import Seq
from cdskit.pad import process_record_padding
from content_digest_cache import cached_sha256_file
from fasta_sequence_store import fasta_records
from gff2genestat import disable_structure, process_single_gff

POLICY_VERSION = 1
GFF_COLUMNS = ['sequence', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
TRAIT_COLUMNS = ['gene_id', 'feature_size', 'num_intron', 'intron_positions', 'chromosome', 'start', 'end',
                 'strand', 'feature_blocks', 'feature_type', 'gff_transcript_id', 'utr_blocks', 'cds_first_phase',
                 'splice_mode', 'feature_block_sequences', 'feature_block_strands', 'transcript_junction_positions',
                 'phase_status', 'cds_partial', 'utr_status', 'structure_status']


def digest(path):
    return cached_sha256_file(Path(path))[0]


def sequence_digest(sequence):
    return hashlib.sha256(sequence.encode('ascii')).hexdigest()


def evaluate(sequence, code=1, phase=None):
    sequence = re.sub(r'\s+', '', sequence).upper()
    if not sequence or re.search(r'[^ACGTRYSWKMBDHVN]', sequence):
        return {'accepted': False, 'reason': 'empty_or_invalid_nucleotide', 'sequence': sequence}
    table = CodonTable.unambiguous_dna_by_id[int(code)]
    stops = set(table.stop_codons) - set(table.forward_table)
    dual = set(table.stop_codons) & set(table.forward_table)
    if phase is None:
        # CDSKit chooses padding. GeneGalleon independently evaluates whether
        # that result meets the scientific admission policy; padding is not
        # evidence that the input is a functional gene.
        padded = str(process_record_padding('candidate', sequence, int(code), 'N')['new_seq'])
        placements = [(h, len(padded) - len(sequence) - h) for h in range(3)
                      if 0 <= len(padded) - len(sequence) - h <= 2 and
                      padded == 'N' * h + sequence + 'N' * (len(padded) - len(sequence) - h)]
        if not placements:
            raise ValueError('CDSKit padding changed sequence content')
        head, tail = placements[0]
    else:
        # GFF phase constrains the first incomplete codon; never trim each
        # exon or ask an unconstrained frame search to override this evidence.
        head = (3 - int(phase)) % 3
        tail = (-(head + len(sequence))) % 3
        padded = 'N' * head + sequence + 'N' * tail
    internal_codons = [padded[i:i + 3] for i in range(0, len(padded) - 3, 3)]
    stop_count = sum(codon in stops for codon in internal_codons)
    reason = 'accepted'
    if stop_count:
        reason = 'premature_stop'
    elif any(re.search('[^ACGT]', codon) for codon in internal_codons[1 if head else 0:]):
        reason = 'internal_ambiguous_bases'
    elif any(codon in dual for codon in internal_codons):
        reason = 'translation_uncertain'
    elif head and phase is None:
        zero_frames = []
        for candidate_head in range(3):
            candidate = 'N' * candidate_head + sequence
            candidate += 'N' * ((-len(candidate)) % 3)
            if not any(candidate[i:i + 3] in stops for i in range(0, len(candidate) - 3, 3)):
                zero_frames.append(candidate_head)
        if len(zero_frames) != 1:
            reason = 'ambiguous_frame'
    return {'accepted': reason == 'accepted', 'reason': reason, 'sequence': padded,
            'source_length': len(sequence), 'head_padding': head, 'tail_padding': tail,
            'frame_offset': (-head) % 3, 'internal_stop_count': stop_count,
            'genetic_code': int(code), 'phase_constraint': phase,
            'source_sha256': sequence_digest(sequence), 'selected_sha256': sequence_digest(padded)}


def extension_related(first, second):
    """Only use length to rank an exact, in-frame extension of the same sequence."""
    short, long = sorted((first, second), key=len)
    if len(short) == len(long):
        return short == long
    start = long.find(short)
    return start >= 0 and start % 3 == 0 and (len(long) - len(short)) % 3 == 0


def select_candidate(supplied, derived=None, code=1, phase=None):
    candidates = {'supplied': evaluate(supplied, code)}
    if derived is not None:
        candidates['gff_genome'] = evaluate(derived, code, phase)
    valid = [name for name, value in candidates.items() if value['accepted']]
    if not valid:
        return None, 'no_valid_cds', candidates
    chosen = valid[0]  # supplied is the default authority when both are valid
    reason = 'supplied_valid' if chosen == 'supplied' else 'gff_genome_rescue'
    if len(valid) == 2:
        a, b = (candidates[name]['sequence'] for name in valid)
        if a == b:
            reason = 'sources_agree'
        elif extension_related(a, b):
            chosen = max(valid, key=lambda name: len(candidates[name]['sequence']))
            reason = 'longer_in_frame_extension'
        else:
            reason = 'supplied_valid_unresolved_source_disagreement'
    return chosen, reason, candidates


def extract_genomic_candidates(traits, genome):
    requests = {}
    pieces = {}
    for row in traits.itertuples(index=False):
        blocks = str(row.feature_blocks).split(';')
        seqids = [unquote(x) for x in str(row.feature_block_sequences).split(';')]
        strands = str(row.feature_block_strands).split(';')
        pieces[row.gene_id] = [None] * len(blocks)
        for index, (block, seqid, strand) in enumerate(zip(blocks, seqids, strands, strict=True)):
            start, end = map(int, block.split('-'))
            requests.setdefault(seqid, []).append((row.gene_id, index, start, end, strand))
    seen = set()
    # Stream contigs, rather than retaining an entire plant genome in memory.
    for identifier, _header, sequence in fasta_records(Path(genome)):
        if identifier not in requests:
            continue
        if identifier in seen:
            raise ValueError(f'Duplicate reference contig: {identifier}')
        seen.add(identifier)
        for gene, index, start, end, strand in requests[identifier]:
            if start < 1 or end > len(sequence):
                raise ValueError(f'CDS outside reference bounds: {gene}')
            part = sequence[start - 1:end].upper()
            pieces[gene][index] = str(Seq(part).reverse_complement()) if strand == '-' else part
    return {gene: ''.join(parts) for gene, parts in pieces.items() if all(part is not None for part in parts)}


def resolve_records(records, traits, genomic, code):
    indices = {row.gene_id: index for index, row in traits.iterrows()}
    output, decisions = [], []
    for identifier, header, sequence in records:
        index = indices.get(identifier)
        phase = None
        if index is not None and traits.at[index, 'phase_status'] == 'consistent':
            phase = int(traits.at[index, 'cds_first_phase'])
        selected, reason, candidates = select_candidate(sequence, genomic.get(identifier), code, phase)
        entry = {'gene_id': identifier, 'selected_source': selected, 'reason': reason, 'candidates': {}}
        for name, candidate in candidates.items():
            entry['candidates'][name] = {k: v for k, v in candidate.items() if k != 'sequence'}
        if selected is None:
            if index is not None:
                disable_structure(traits, index, 'excluded_cds')
            entry['structure_status'] = 'excluded_cds'
        else:
            chosen = candidates[selected]
            output.append((identifier, header, chosen['sequence']))
            entry['selected_sha256'] = chosen['selected_sha256']
            entry['structure_status'] = 'unverified_source'
            if index is not None:
                raw = genomic.get(identifier)
                # Any 5' padding changes intron offsets. Keep the sequence but
                # do not silently reuse the old coordinates or coding phases.
                compatible = raw is not None and chosen['sequence'] == raw + 'N' * ((-len(raw)) % 3)
                if compatible:
                    traits.at[index, 'structure_status'] = 'sequence_verified'
                    entry['structure_status'] = 'sequence_verified'
                else:
                    disable_structure(traits, index, 'sequence_not_coordinate_matched')
                    entry['structure_status'] = 'sequence_not_coordinate_matched'
        decisions.append(entry)
    return output, decisions


def resolve(cds, gff, genome, output_dir, code=1):
    cds, output_dir = Path(cds), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / cds.name
    report_path = output_dir / (cds.name + '.resolution.json')
    traits_path = output_dir / (cds.name + '.traits.tsv')
    meta_path = output_dir / (cds.name + '.resolution.meta.json')
    if any(destination.resolve() == Path(path).resolve() for path in (cds, gff, genome) if path):
        raise ValueError('Resolved output must not overwrite a source file')
    sources = {name: {'path': str(Path(path).resolve()), 'sha256': digest(path)}
               for name, path in [('cds', cds), ('gff', gff), ('genome', genome)] if path}
    support = Path(__file__).parent
    contract = {'policy_version': POLICY_VERSION, 'genetic_code': int(code), 'sources': sources,
                'dependencies': {name: importlib.metadata.version(name) for name in ['cdskit', 'biopython']},
                'implementation': {name: digest(support / name) for name in
                                   ['cds_resolution.py', 'gff2genestat.py', 'gff_feature_structure.py']}}
    with (output_dir / (cds.name + '.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if meta_path.exists():
            old = json.loads(meta_path.read_text())
            if old.get('contract') == contract and destination.exists() and traits_path.exists() and report_path.exists():
                if old.get('output_sha256') == digest(destination) and old.get('traits_sha256') == digest(traits_path):
                    return destination, traits_path, report_path
        records = list(fasta_records(cds))
        if len({record[0] for record in records}) != len(records):
            raise ValueError('Duplicate CDS identifier')
        traits = pd.DataFrame(columns=TRAIT_COLUMNS)
        if gff:
            traits = process_single_gff(Path(gff).name, str(Path(gff).parent), [r[0] for r in records],
                                        'CDS', 'longest', GFF_COLUMNS, TRAIT_COLUMNS, 'report', 'report')
            if traits.empty and records:
                raise ValueError('No GFF IDs match supplied CDS; cannot resolve source identity')
        genomic = extract_genomic_candidates(traits, genome) if genome and not traits.empty else {}
        resolved, decisions = resolve_records(records, traits, genomic, code)
        if any(digest(entry['path']) != entry['sha256'] for entry in sources.values()):
            raise ValueError('CDS resolution source changed during validation')
        # A nonempty source with no usable genes is explicit failure, not an
        # apparently successful empty species consumed by downstream tools.
        if records and not resolved:
            raise ValueError('No usable CDS remains after translation validation')
        with tempfile.TemporaryDirectory(prefix='.cds-resolution-', dir=output_dir) as temporary:
            staged = Path(temporary) / destination.name
            opener = gzip.open if destination.name.endswith('.gz') else open
            with opener(staged, 'wt') as handle:
                for identifier, _header, sequence in resolved:
                    handle.write(f'>{identifier}\n{sequence}\n')
            selected_hashes = {identifier: sequence_digest(sequence) for identifier, _, sequence in resolved}
            traits = traits.loc[traits['gene_id'].isin(selected_hashes)].copy()
            traits['cds_sequence_sha256'] = traits['gene_id'].map(selected_hashes)
            staged_traits = Path(temporary) / 'traits.tsv'
            traits.to_csv(staged_traits, sep='\t', index=False)
            report = {'contract': contract, 'output_sha256': digest(staged), 'traits_sha256': digest(staged_traits),
                      'input_count': len(records), 'accepted_count': len(resolved),
                      'excluded_count': len(records) - len(resolved), 'decisions': decisions}
            staged_report = Path(temporary) / 'report.json'
            staged_report.write_text(json.dumps(report, sort_keys=True) + '\n')
            os.replace(staged, destination)
            os.replace(staged_traits, traits_path)
            os.replace(staged_report, report_path)
            staged_meta = Path(temporary) / 'meta.json'
            staged_meta.write_text(json.dumps({k: v for k, v in report.items() if k != 'decisions'}, sort_keys=True) + '\n')
            os.replace(staged_meta, meta_path)  # commit marker published last
    return destination, traits_path, report_path


def resolved_view(source_dir, resolution_dir, view_dir):
    """Bind downstream inputs to verified resolutions; preserve unprocessed species."""
    source_dir, resolution_dir, view_dir = map(Path, (source_dir, resolution_dir, view_dir))
    view_dir.mkdir(parents=True, exist_ok=True)
    selected = {}
    for source in sorted(source_dir.iterdir()):
        if not source.is_file() or not re.search(r'\.(fa|fas|fasta|fna)(\.gz)?$', source.name):
            continue
        target = source.resolve()
        report_path = resolution_dir / (source.name + '.resolution.meta.json')
        if report_path.exists():
            report = json.loads(report_path.read_text())
            for entry in report['contract']['sources'].values():
                if digest(entry['path']) != entry['sha256']:
                    raise ValueError(f'Stale CDS resolution; rerun genome annotation: {source.name}')
            if report['contract']['sources']['cds']['path'] != str(source.resolve()):
                raise ValueError(f'CDS resolution source changed: {source.name}')
            target = (resolution_dir / source.name).resolve()
            if digest(target) != report['output_sha256']:
                raise ValueError(f'Changed resolved CDS: {source.name}')
        selected[source.name] = {'path': str(target), 'sha256': digest(target)}
    identity = hashlib.sha256(json.dumps(selected, sort_keys=True).encode()).hexdigest()
    destination = view_dir / identity
    if destination.exists():
        for name, entry in selected.items():
            if digest(destination / name) != entry['sha256']:
                raise ValueError(f'Changed CDS input view: {name}')
        return destination
    # Regular hardlinks keep the existing FASTA-store/no-symlink contract.
    # Each view has an immutable content identity; no user source is rewritten.
    with tempfile.TemporaryDirectory(prefix='.view-', dir=view_dir) as temporary:
        staged = Path(temporary) / 'files'
        staged.mkdir()
        for name, entry in selected.items():
            os.link(entry['path'], staged / name)
        try:
            os.rename(staged, destination)
        except OSError:
            if not destination.is_dir():
                raise
            for name, entry in selected.items():
                if digest(destination / name) != entry['sha256']:
                    raise ValueError(f'Concurrent CDS view differs: {name}') from None
    return destination


def backup_family_outputs(workspace, family):
    """Preserve manifest-bound derived files before an explicit policy rebuild."""
    workspace = Path(workspace).resolve()
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]+', family):
        raise ValueError('Invalid family identity')
    inventory = {}
    for manifest in sorted((workspace / 'output/artifact_provenance/genome_annotation').glob(family + '.*.json')):
        payload = json.loads(manifest.read_text())
        if payload.get('family_id') != family:
            raise ValueError('Backup manifest family changed')
        inventory[str(manifest.relative_to(workspace))] = digest(manifest)
        for entry in payload.get('outputs', []) + payload.get('optional_outputs', []):
            if entry.get('scope') != 'workspace' or entry.get('artifact_type') != 'file':
                raise ValueError('Unrecognized output in CDS policy migration')
            path = workspace / entry['path']
            if not path.exists():
                continue
            if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(workspace / 'output'):
                raise ValueError('Unsafe output in CDS policy migration')
            inventory[str(path.relative_to(workspace))] = digest(path)
    if not inventory:
        return None
    identity = hashlib.sha256(json.dumps(inventory, sort_keys=True).encode()).hexdigest()
    parent = workspace / 'output/cds_resolution_history' / family
    parent.mkdir(parents=True, exist_ok=True)
    destination = parent / identity
    if destination.exists():
        for relative, expected in inventory.items():
            if digest(destination / relative) != expected:
                raise ValueError('Changed CDS policy recovery copy')
        return destination
    with tempfile.TemporaryDirectory(prefix='.backup-', dir=parent) as temporary:
        staged = Path(temporary) / 'files'
        staged.mkdir()
        for relative, expected in inventory.items():
            target = staged / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(workspace / relative, target)
            if digest(target) != expected or digest(workspace / relative) != expected:
                raise ValueError('Output changed during CDS policy backup')
        (staged / 'inventory.json').write_text(json.dumps(inventory, sort_keys=True) + '\n')
        os.rename(staged, destination)
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cds')
    parser.add_argument('--gff')
    parser.add_argument('--genome')
    parser.add_argument('--output-dir')
    parser.add_argument('--backup-family')
    parser.add_argument('--workspace')
    parser.add_argument('--genetic-code', type=int, default=1)
    parser.add_argument('--source-dir')
    parser.add_argument('--view-dir')
    args = parser.parse_args()
    if args.backup_family:
        if not args.workspace:
            parser.error('--backup-family requires --workspace')
        backup_family_outputs(args.workspace, args.backup_family)
        return
    if not args.output_dir:
        parser.error('--output-dir is required')
    if args.source_dir:
        if not args.view_dir:
            parser.error('--source-dir requires --view-dir')
        print(resolved_view(args.source_dir, args.output_dir, args.view_dir))
    else:
        if not args.cds:
            parser.error('--cds is required')
        resolve(args.cds, args.gff, args.genome, args.output_dir, args.genetic_code)


if __name__ == '__main__':
    main()

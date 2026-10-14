#!/usr/bin/env python3
"""Annotation-masking admission comparison using independent GFF and donors.

The fixed native ANNEVO prediction must precede masking and be bound to the same
genome. Donor miniprot placement sees only DNA/proteins, never held-out labels.
Reference absence is not a biological negative; off-reference additions are
reported separately. This measures conservative model admission, not a complete
cohort analysis or calibrated annotation accuracy.
"""
import argparse
import copy
import hashlib
import json
import random
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
import annevo_gene_models as annevo
import gene_model_refinement as refinement
import rescue_gene_models as rescue
from fasta_sequence_store import fasta_records, open_text
from format_species_annotation.common import parse_gff_attributes
from gene_model_catalog import build_catalog, indexed_genome


def attributes(text):
    return {key: values for key, values in parse_gff_attributes(text).items()}


def reference_rows(path, genome):
    """Keep complete original gene envelopes in the exact tested contigs."""
    rows, graph, genes = [], {}, set()
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\r\n').split('\t')
            if len(f) != 9:
                raise ValueError('Malformed independent reference GFF')
            if f[0] not in genome.references:
                continue
            attrs = attributes(f[8])
            identifier = attrs.get('ID', (None,))[0]
            parents = attrs.get('Parent', ())
            if identifier:
                previous = graph.get(identifier)
                if previous and previous != (f[2], parents):
                    raise ValueError('Conflicting independent reference feature identity')
                graph[identifier] = (f[2], parents)
                if f[2] == 'gene' and 0 <= int(f[3]) - 1 < int(f[4]) <= genome.get_reference_length(f[0]):
                    genes.add(identifier)
            rows.append((f, attrs))
    def owners(identifier, seen=()):
        if identifier in seen:
            raise ValueError('Cycle in independent reference GFF')
        if identifier in genes:
            return {identifier}
        return set().union(*(owners(p, (*seen, identifier)) for p in graph.get(identifier, ('', ()))[1]))
    kept = []
    for f, attrs in rows:
        roots = set().union(*(owners(i) for i in attrs.get('ID', ()))) if attrs.get('ID') else set()
        roots |= set().union(*(owners(p) for p in attrs.get('Parent', ())))
        if roots and roots <= genes:
            if len(roots) != 1:
                raise ValueError('Shared reference feature has ambiguous gene ownership')
            kept.append((f, attrs, next(iter(roots))))
    return kept


def write_annotation(directory, rows, genome, *, hidden_genes=(), hidden_blocks=None):
    directory.mkdir(parents=True, exist_ok=False)
    hidden_blocks = hidden_blocks or {}
    kept, coding = [], defaultdict(list)
    for f, attrs, owner in rows:
        block = int(f[3]) - 1, int(f[4])
        if owner in hidden_genes or (f[2] in {'CDS', 'exon'} and block in hidden_blocks.get(owner, set())):
            continue
        kept.append((f, attrs, owner))
        if f[2] == 'CDS':
            for parent in attrs.get('Parent', ()):
                coding[parent].append(f)
    gff, cds = directory / 'reference.gff3', directory / 'reference.cds.fa'
    gff.write_text('##gff-version 3\n' + ''.join('\t'.join(f) + '\n' for f, _, _ in kept))
    with cds.open('w') as output:
        for transcript, features in sorted(coding.items()):
            strand, seqid = features[0][6], features[0][0]
            blocks = sorted([(int(f[3]) - 1, int(f[4])) for f in features], reverse=strand == '-')
            sequence = ''.join(str(Seq(genome.fetch(seqid, a, b)).reverse_complement()) if strand == '-'
                               else genome.fetch(seqid, a, b) for a, b in blocks)
            output.write('>' + transcript + '\n' + sequence + '\n')
    return gff, cds


def donor_placements(args, directory, genome):
    proteins, queries, removed = {}, {}, Counter()
    protein_file = directory / 'donors.pep.fa'
    with protein_file.open('w') as output:
        for entry in args.donor_cds:
            name, path = entry.split('=', 1)
            if name == args.species or name in proteins:
                raise ValueError('Donor species must be distinct from the target and each other')
            proteins[name] = {}
            for identifier, _, dna in fasta_records(Path(path)):
                dna = dna.upper()
                if len(dna) % 3 or set(dna) - set('ACGT'):
                    removed['frame_or_ambiguity'] += 1
                    continue
                protein = str(Seq(dna).translate()).removesuffix('*')
                if not protein or '*' in protein:
                    removed['internal_stop'] += 1
                    continue
                if identifier in proteins[name]:
                    raise ValueError('Duplicate donor CDS identifier')
                proteins[name][identifier] = protein
                query = 'q' + hashlib.sha256((name + '\0' + identifier).encode()).hexdigest()[:24]
                if query in queries:
                    raise ValueError('Donor query identifier collision')
                queries[query] = (name, identifier)
                output.write('>' + query + '\n' + protein + '\n')
    if len(proteins) < 2:
        raise ValueError('At least two genuinely independent donor species are required')
    index = directory / 'target.mpi'
    raw_gff = directory / 'donor_placements.gff3'
    commands = [['miniprot', '-T', '1', '-t', str(args.cpus), '-d', str(index), str(args.genome)],
                ['miniprot', '-u', '-t', str(args.cpus), '-G', '20000', '-N30', '--outs=0.5', '--gff', str(index), str(protein_file)]]
    for number, command in enumerate(commands):
        with (directory / ('miniprot_' + str(number) + '.stderr')).open('w') as err:
            if number:
                with raw_gff.open('w') as output:
                    subprocess.run(command, stdout=output, stderr=err, check=True)
            else:
                subprocess.run(command, stderr=err, check=True)
    models, regions = [], []
    for number, raw in enumerate(rescue.read_miniprot(raw_gff)):
        donor, query = queries[raw['query']]
        a, b = min(x[0] for x in raw['cds']), max(x[1] for x in raw['cds'])
        region = dict(id='placement' + str(number), donor=donor, query=query, target=args.species,
                      seqid=raw['seqid'], start=max(0, a - 300), end=min(genome.get_reference_length(raw['seqid']), b + 300),
                      genome_only=True, annevo_bounded_placement=True)
        raw.update(query=region['id'], evidence=region, search='genome_fallback')
        # Match the production worker: complete real donor-supported partial
        # termini before applying placement/ownership flags. The previous
        # assay skipped this stage, understating whole-path confirmations.
        checked = rescue.validate_model(raw, genome, 1, refinement.DEFAULTS)
        checked = rescue.complete_terminals(checked, genome, 1, refinement.DEFAULTS,
                                            proteins[donor][query], rescue.validate_model)
        checked = rescue.check_interval(checked)
        models.append(checked)
        if b - a <= refinement.DEFAULTS['max_interval']:
            regions.append(region)
    return models, regions, proteins, {'commands': commands, 'donor_queries': len(queries),
                                       'donor_species': sorted(proteins), 'excluded_donor_CDS': dict(removed)}


def shape(model):
    blocks = model['blocks'] if 'blocks' in model else model['cds']
    return model['seqid'], model['strand'], tuple(tuple(b) for b in blocks)


def sampled_masks(truth, exon_choices, fraction, seed):
    """Freeze label masks independently of process hash/set iteration order."""
    rng = random.Random(seed)
    hidden = set(rng.sample(sorted(truth), max(1, int(len(truth) * fraction))))
    exons = set(rng.sample(sorted(exon_choices), max(1, int(len(exon_choices) * fraction)))) if exon_choices else set()
    return hidden, {g: {rng.choice(exon_choices[g])} for g in sorted(exons)}


def evaluate_rescue(models, gff, truth, hidden, species, audit_path=None):
    with gff.open() as stream:
        existing = rescue.original_ownership(stream, rescue.attributes)
    ownership = rescue.OwnershipIndex(existing)
    records = copy.deepcopy(models)
    rescue.reassess_unanchored_models(records, 2, ownership=ownership)
    records = rescue.consolidate(records, existing, species)
    if audit_path:
        rescue.atomic_json(audit_path, records)
    accepted = [m for m in records if m['status'] == 'accepted']
    truth_shapes = {shape(c): gene for gene, c in truth.items()}
    recovered = {truth_shapes[shape(m)] for m in accepted if shape(m) in truth_shapes and truth_shapes[shape(m)] in hidden}
    nonexact_in_hidden = [m for m in accepted if shape(m) not in truth_shapes
                         and any(m['seqid'] == truth[g]['seqid'] and m['strand'] == truth[g]['strand']
                                 and any(a[0] < b[1] and a[1] > b[0] for a in m['cds'] for b in truth[g]['blocks']) for g in hidden)]
    return {'hidden_genes': len(hidden), 'exact_recovered_genes': len(recovered),
            'nonexact_additions_at_hidden_loci': len(nonexact_in_hidden), 'accepted_additions': len(accepted),
            'off_reference_additions_unknown': sum(not any(m['seqid'] == c['seqid'] and m['strand'] == c['strand']
                                                         and any(a[0] < b[1] and a[1] > b[0] for a in m['cds'] for b in c['blocks'])
                                                         for c in truth.values()) for m in accepted),
            'ANNEVO_supported_additions': sum(bool(m.get('evidence', {}).get('annevo') or any(s.get('annevo') for s in m.get('support', []))) for m in accepted)}


def evaluate_refinement(models, catalog, truth, hidden, species, audit_path=None):
    loci = {g['gene_id']: g for g in catalog['loci']}
    annotated = []
    edges, seen, query_loci = [], set(), defaultdict(set)
    for raw in models:
        owned = [g for g in loci.values() if g['seqid'] == raw['seqid'] and g['strand'] == raw['strand']
                 and any(a[0] < b[1] and a[1] > b[0] for c in g['candidates'] for a in c['blocks'] for b in raw['cds'])]
        if len(owned) != 1:
            continue
        gene = owned[0]['gene_id']
        evidence = raw['evidence']
        model = copy.deepcopy(raw)
        model.update(gene_id=gene, donor_species=evidence['donor'], donor_candidate=evidence['query'])
        model['problems'] = [p for p in model['problems'] if p != 'unanchored_genome_search']
        if evidence.get('annevo'):
            # Refinement has an existing, uniquely owned annotated locus. Its
            # confirmation scope differs from unanchored new-gene rescue.
            model['evidence']['genome_only'] = False
        annotated.append(model)
        key = gene, evidence['donor'], evidence['query']
        if key not in seen:
            seen.add(key)
            edges.append(dict(species_a=species, gene_a=gene, species_b=evidence['donor'], gene_b=evidence['query'], ambiguous=False))
    annevo.confirm_transferred_structure([m for m in annotated if m.get('evidence', {}).get('annevo')],
                                         [m for m in annotated if not m.get('evidence', {}).get('annevo')])
    for model in annotated:
        if not model['problems']:
            query_loci[model['donor_species'], model['donor_candidate']].add(model['gene_id'])
    for edge in edges:
        edge['ambiguous'] = len(query_loci[edge['species_b'], edge['gene_b']]) > 1
    predictions = refinement.classify_predictions(annotated, catalog, edges, refinement.DEFAULTS, [], catalog['sources']['genome']['sha256'])
    if audit_path:
        rescue.atomic_json(audit_path, predictions)
    accepted = [r for r in predictions if r['status'] == 'accepted']
    recovered = {r['gene_id'] for r in accepted if r['gene_id'] in hidden
                 and shape(r['candidate']) == shape(truth[r['gene_id']])}
    intact = [r for r in accepted if r['gene_id'] in truth and r['gene_id'] not in hidden
              and r['candidate']['quality'].get('representative_eligible')]
    return {'hidden_exon_genes': len(hidden), 'exact_repair_admissions': len(recovered),
            'nonexact_repair_admissions': sum(r['gene_id'] in hidden and shape(r['candidate']) != shape(truth[r['gene_id']]) for r in accepted),
            'exact_representative_eligible_repairs': len({r['gene_id'] for r in accepted if r['gene_id'] in hidden
                                                        and shape(r['candidate']) == shape(truth[r['gene_id']])
                                                        and r['candidate']['quality'].get('representative_eligible')}),
            'intact_locus_representative_admissions': len(intact),
            'ANNEVO_supported_admissions': sum(bool(r['candidate'].get('support', {}).get('annevo')) for r in accepted)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('genome', 'reference-gff', 'annevo-gff', 'receipt', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--species', default='Aethionema_arabicum')
    parser.add_argument('--donor-cds', action='append', required=True, help='SPECIES=CDS_FASTA; real distinct donor species')
    parser.add_argument('--seed', action='append', type=int)
    parser.add_argument('--fraction', type=float, default=0.25)
    parser.add_argument('--cpus', type=int, default=4)
    parser.add_argument('--mask-manifest', type=Path, help='Replay source-bound gene/exon masks for equivalent comparisons')
    args = parser.parse_args()
    if not 0 < args.fraction < 1 or args.cpus < 1:
        raise ValueError('Invalid holdout/resource bound')
    args.output.mkdir(parents=True, exist_ok=False)
    input_paths = [args.genome, args.reference_gff, args.annevo_gff, args.receipt,
                   *(Path(v.split('=', 1)[1]) for v in args.donor_cds)]
    hashes = rescue.digest_paths(input_paths)
    replay = None
    if args.mask_manifest:
        manifest_sha = rescue.digest(args.mask_manifest)
        replay = json.loads(args.mask_manifest.read_text())
        if rescue.digest(args.mask_manifest) != manifest_sha:
            raise ValueError('Mask manifest changed during reading')
        if (replay.get('schema') != 1 or replay.get('genome_sha256') != hashes[str(args.genome)]
                or replay.get('reference_gff_sha256') != hashes[str(args.reference_gff)]):
            raise ValueError('Mask manifest differs from the benchmark genome/reference')
        hashes[str(args.mask_manifest)] = manifest_sha
    results = []
    with indexed_genome(args.genome) as genome:
        rows = reference_rows(args.reference_gff, genome)
        ref_gff, ref_cds = write_annotation(args.output / 'complete_reference', rows, genome)
        full = build_catalog(args.species, ref_cds, ref_gff, args.genome)
        truth = {g['gene_id']: min((c for c in g['candidates'] if c['quality']['valid_orf']), key=lambda c: c['candidate_id'])
                 for g in full['loci'] if any(c['quality']['valid_orf'] for c in g['candidates'])}
        source_ids = {g['gene_id']: g['candidates'][0]['source_gene_id'] for g in full['loci']}
        table = args.output / 'native_candidates.tsv'
        table.write_text('species\tgff\treceipt\n' + '\t'.join((args.species, str(args.annevo_gff.resolve()), str(args.receipt.resolve()))) + '\n')
        frozen = annevo.freeze_inputs(table, {args.species: {'genome': str(args.genome), 'genetic_code': 1}})
        baseline, regions, proteins, placement = donor_placements(args, args.output, genome)
        native_models, native_audit = annevo.nominate(frozen, args.species, regions, proteins, genome, refinement.DEFAULTS)
        checked = []
        for model in native_models:
            m = rescue.validate_model(model, genome, 1, refinement.DEFAULTS)
            m['problems'] = sorted(set(m['problems']) | set(model['annevo_problems']))
            checked.append(rescue.check_interval(m))
        annevo.confirm_transferred_structure(checked, baseline)
        rescue.atomic_json(args.output / 'annevo_nomination.json', native_audit)
        exon_choices = {g: [tuple(b[:2]) for b in c['blocks'][1:-1] if (b[1] - b[0]) % 3]
                        for g, c in truth.items() if len(c['blocks']) >= 3}
        exon_choices = {g: v for g, v in exon_choices.items() if v}
        masks = []
        seeds = args.seed or ([r['seed'] for r in replay['seeds']] if replay else [20261014, 20261015, 20261016])
        if len(set(seeds)) != len(seeds):
            raise ValueError('Duplicate mask seed')
        for seed in seeds:
            if replay:
                matching = [r for r in replay['seeds'] if r['seed'] == seed]
                if len(matching) != 1:
                    raise ValueError('Missing or duplicate replay seed')
                selected = matching[0]
                hidden = set(selected['hidden_genes'])
                hidden_exons = set(selected['hidden_exon_blocks'])
                if (not hidden or not hidden <= truth.keys() or not hidden_exons <= exon_choices.keys()
                        or len(hidden) != len(selected['hidden_genes'])):
                    raise ValueError('Unknown or duplicate replay locus')
                masked = {}
                for g, blocks in selected['hidden_exon_blocks'].items():
                    if len(blocks) != 1 or tuple(blocks[0]) not in exon_choices[g]:
                        raise ValueError('Replay exon is not an eligible original internal coding block')
                    masked[source_ids[g]] = {tuple(blocks[0])}
            else:
                hidden, sampled = sampled_masks(truth, exon_choices, args.fraction, seed)
                hidden_exons = set(sampled)
                masked = {source_ids[g]: blocks for g, blocks in sampled.items()}
            masks.append({'seed': seed, 'hidden_genes': sorted(hidden),
                          'hidden_exon_blocks': {g: [list(b) for b in sorted(masked[source_ids[g]])] for g in sorted(hidden_exons)}})
            gff, _ = write_annotation(args.output / ('gene_holdout_' + str(seed)), rows, genome,
                                       hidden_genes={source_ids[g] for g in hidden})
            gene = {'seed': seed, 'without_ANNEVO': evaluate_rescue(baseline, gff, truth, hidden, args.species,
                                                                  gff.parent / 'without_ANNEVO_admission.json'),
                    'with_ANNEVO': evaluate_rescue(baseline + checked, gff, truth, hidden, args.species,
                                                   gff.parent / 'with_ANNEVO_admission.json')}
            gff, cds = write_annotation(args.output / ('exon_holdout_' + str(seed)), rows, genome, hidden_blocks=masked)
            catalog = build_catalog(args.species, cds, gff, args.genome)
            catalog['annotation_spans'] = refinement.annotation_ownership_spans(gff, catalog)
            exon = {'seed': seed, 'without_ANNEVO': evaluate_refinement(baseline, catalog, truth, hidden_exons, args.species,
                                                                       gff.parent / 'without_ANNEVO_admission.json'),
                    'with_ANNEVO': evaluate_refinement(baseline + checked, catalog, truth, hidden_exons, args.species,
                                                        gff.parent / 'with_ANNEVO_admission.json')}
            results.append({'gene_masking': gene, 'exon_masking': exon})
    if rescue.digest_paths(hashes) != hashes:
        raise ValueError('Benchmark source inputs changed')
    rescue.atomic_json(args.output / 'mask_manifest.json', {'schema': 1, 'genome_sha256': hashes[str(args.genome)],
                      'reference_gff_sha256': hashes[str(args.reference_gff)], 'seeds': masks})
    report = {'scope': 'one genome interval; conservative admission, not final cohort/orthology or biological truth',
              'independent_reference': True, 'synthetic_RNA_evidence': False, 'reference_genes': len(full['loci']),
              'intact_reference_genes': len(truth), 'exon_maskable_genes': len(exon_choices),
              'input_sha256': hashes, 'parameters': refinement.DEFAULTS, 'placement': placement, 'seeds': results,
              'implementation': refinement.implementation(), 'dependencies': refinement.dependency_identities(),
              'benchmark_script_sha256': rescue.digest(__file__), 'mask_manifest_sha256': rescue.digest(args.output / 'mask_manifest.json'),
              'miniprot': {'version': subprocess.check_output(['miniprot', '--version'], text=True).strip(),
                          'binary_sha256': rescue.digest(shutil.which('miniprot'))},
              'ANNEVO_prediction_reused_without_rerunning_GPU': True,
              'annotation_absence_is_not_false_positive_truth': True}
    rescue.atomic_json(args.output / 'comparison.json', report)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

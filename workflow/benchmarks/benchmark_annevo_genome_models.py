#!/usr/bin/env python3
"""Whole-genome, source-bound masking and negative-control admission audit.

This uses actual miniprot placement, donor-supported terminal completion,
ANNEVO genomic confirmation, ownership and classification. Correspondence is
an explicitly approximate placement graph, not independently verified synteny.
All reference isoforms count as matches; absence is unknown, not error truth.
The native-experiment option validates original external evidence without
manufacturing a GeneGalleon native receipt or relaxing production import rules.
"""
import argparse
import copy
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import benchmark_annevo_masked_models as bm
from gene_model_catalog import build_catalog, indexed_genome
from gene_model_selection import select_representatives

ann, ref, rescue = bm.annevo, bm.refinement, bm.rescue


def external_prediction(directory, species, genome, native_gff):
    directory = Path(directory)
    paths = {n: directory / n for n in ('manifest.json', 'upstream_exit.json', 'runtime_inventory.json')}
    hashes = rescue.digest_paths(paths.values())
    manifest, exit_record, inventory = (json.loads(paths[n].read_text()) for n in paths)
    inputs, runtime = manifest['inputs'], manifest['runtime']
    if (type(manifest.get('schema_version')) is not int or manifest.get('schema_version') != 1 or manifest.get('kind') != 'predictor-native-package-v1'
            or manifest.get('tool') != 'ANNEVO' or manifest.get('species') != species
            or manifest['options'].get('lineage') != 'Magnoliopsida'
            or manifest['files'].get(inputs['genome']) != rescue.digest(genome)
            or manifest['files'].get(runtime['inventory']) != hashes[str(paths['runtime_inventory.json'])]
            or inventory.get('source_commit') != runtime.get('source_commit')
            or not re.fullmatch('[0-9a-f]{40}', str(runtime.get('source_commit')))
            or exit_record.get('manifest_sha256') != hashes[str(paths['manifest.json'])]
            or type(exit_record.get('upstream_returncode')) is not int or exit_record['upstream_returncode'] != 0
            or type(exit_record.get('integration_returncode')) is not int or exit_record['integration_returncode'] != 0
            or exit_record.get('gpu_probe', {}).get('gpu_verified') is not True
            or exit_record['gpu_probe'].get('tool') != 'ANNEVO'
            or exit_record.get('artifact_sha256') != rescue.digest(native_gff)
            or type(exit_record.get('artifact_bytes')) is not int
            or exit_record.get('artifact_bytes') != Path(native_gff).stat().st_size
            or exit_record.get('upstream_completion_reported') is not True):
        raise ValueError('External native experiment failed its source/GPU/artifact binding')
    command = exit_record['command']
    for option, expected in [('-g', inputs['genome']), ('-m', runtime['model']),
                             ('-l', 'Magnoliopsida'), ('-o', manifest['root'] + '/prediction.gff3')]:
        if command.count(option) != 1 or command[command.index(option) + 1] != expected:
            raise ValueError('External native command does not match its sealed inputs/model/output')
    model_path = str(Path(runtime['model']).relative_to(inventory['root']))
    models = [row for row in inventory['files'] if row['path'] == model_path]
    if len(models) != 1 or not re.fullmatch('[0-9a-f]{64}', str(models[0].get('sha256'))):
        raise ValueError('External native experiment lacks its sealed model')
    if rescue.digest_paths(hashes) != hashes:
        raise ValueError('External native evidence changed during reading')
    identity = dict(provenance_kind='verified_external_native_experiment',
                    manifest_sha256=hashes[str(paths['manifest.json'])],
                    receipt_sha256=hashes[str(paths['upstream_exit.json'])],
                    runtime_digest=hashes[str(paths['runtime_inventory.json'])],
                    source_commit=runtime['source_commit'], model_sha256=models[0]['sha256'],
                    model_profile='Magnoliopsida', genome_sha256=rescue.digest(genome),
                    GFF_sha256=rescue.digest(native_gff), GFF_bytes=Path(native_gff).stat().st_size)
    # This benchmark-only frozen nomination object retains the original evidence.
    # It is deliberately not accepted by the production freeze_inputs API.
    return {'table': str(paths['manifest.json']), 'files': {**hashes, str(native_gff): identity['GFF_sha256']},
            'predictions': {species: {'gff': str(native_gff), 'receipt': str(paths['upstream_exit.json']), 'identity': identity}}}


def reference_index(catalog):
    refs = {g['gene_id']: g['candidates'] for g in catalog['loci']}
    features = []
    shapes = defaultdict(set)
    for gene, candidates in refs.items():
        for c in candidates:
            shapes[bm.shape(c)].add(gene)
        blocks = [b for c in candidates for b in c['blocks']]
        if blocks:
            c = candidates[0]
            features.append(dict(gene_id=gene, seqid=c['seqid'], strand=c['strand'], cds=blocks,
                                 start=min(b[0] for b in blocks), end=max(b[1] for b in blocks)))
    return refs, shapes, rescue.OwnershipIndex(features)


def model_bounds(m):
    blocks = m.get('blocks', m.get('cds'))
    return {**m, 'cds': blocks, 'start': min(b[0] for b in blocks), 'end': max(b[1] for b in blocks)}


def match_category(model, refs, shapes, index):
    exact = shapes.get(bm.shape(model), set())
    owners = {r['gene_id'] for r in index.overlapping(model_bounds(model))}
    if exact:
        return 'exact_any_reference_isoform', sorted(exact), sorted(owners)
    if len(owners) > 1:
        return 'multiple_reference_loci', [], sorted(owners)
    if not owners:
        return 'off_reference_unknown', [], []
    gene = next(iter(owners))
    blocks = sorted(model.get('blocks', model.get('cds')))
    junctions = tuple((a[1], b[0]) for a, b in zip(blocks, blocks[1:], strict=False))
    if any(junctions == tuple((a[1], b[0]) for a, b in zip(sorted(c['blocks']), sorted(c['blocks'])[1:], strict=False))
           for c in refs[gene]):
        return 'terminal_change_only', [], [gene]
    return 'internal_structure_difference', [], [gene]


def dump_compressed(path, value):
    with gzip.open(path, 'wt', compresslevel=2) as stream:
        json.dump(value, stream, sort_keys=True)


def gene_assay(models, gff, species, hidden, reference, output):
    with gff.open() as stream:
        existing = rescue.original_ownership(stream, rescue.attributes)
    records = copy.deepcopy(models)
    rescue.reassess_unanchored_models(records, 2, ownership=rescue.OwnershipIndex(existing))
    records = rescue.consolidate(records, existing, species)
    accepted = [m for m in records if m['status'] == 'accepted']
    categories, recovered, copies, cases = Counter(), set(), Counter(), []
    for m in accepted:
        category, exact, owners = match_category(m, *reference)
        categories[category] += 1
        recovered.update(set(exact) & hidden)
        copies.update(owners)
        cases.append(dict(model_id=m['model_id'], category=category, exact_reference_loci=exact,
                          overlapping_reference_loci=owners, cds=m['cds'], support=m.get('support', [])))
    dump_compressed(output / 'models.json.gz', records)
    rescue.atomic_json(output / 'accepted_cases.json', cases)
    return dict(hidden_genes=len(hidden), exact_recovered_genes=len(recovered), accepted_additions=len(accepted),
                categories=dict(categories), multiple_models_at_one_reference_locus=sum(n > 1 for n in copies.values()),
                proposal_reasons=dict(Counter(p for m in records if m['status'] != 'accepted' for p in m['problems'])))


def query_locus_map(catalogs):
    aliases = {}
    for catalog in catalogs:
        for g in catalog['loci']:
            for c in g['candidates']:
                for source in c.get('source_cds', []):
                    key = catalog['species'], source['source_fasta_id']
                    value = g['gene_id'], c['candidate_id']
                    if key in aliases and aliases[key] != value:
                        raise ValueError('Ambiguous donor CDS ownership')
                    aliases[key] = value
    return aliases


def refinement_assay(models, catalog, donors, aliases, hidden, reference, output):
    loci = {g['gene_id']: g for g in catalog['loci']}
    _refs, _shapes, ownership = reference_index(catalog)
    annotated, edges, seen = [], [], set()
    query_loci = defaultdict(set)
    for raw in models:
        owned = ownership.overlapping(model_bounds(raw))
        if len(owned) != 1:
            continue
        e = raw['evidence']
        lookup = aliases.get((e['donor'], e['query']))
        if not lookup:
            continue
        gene = owned[0]['gene_id']
        m = copy.deepcopy(raw)
        m.update(gene_id=gene, donor_species=e['donor'], donor_candidate=lookup[1])
        m['evidence']['genome_only'] = False
        m['problems'] = [p for p in m['problems'] if p != 'unanchored_genome_search']
        annotated.append(m)
        if not m['problems']:
            query_loci[e['donor'], lookup[0]].add(gene)
        key = gene, e['donor'], lookup[0]
        if key not in seen:
            seen.add(key)
            edges.append(dict(species_a=catalog['species'], gene_a=gene, species_b=e['donor'], gene_b=lookup[0],
                              ambiguous=False, evidence='approximate_genomic_placement_not_verified_synteny'))
    ann.confirm_transferred_structure([m for m in annotated if m['evidence'].get('annevo')],
                                     [m for m in annotated if not m['evidence'].get('annevo')])
    for edge in edges:
        edge['ambiguous'] = len(query_loci[edge['species_b'], edge['gene_b']]) > 1
    predictions = ref.classify_predictions(annotated, catalog, edges, ref.DEFAULTS, [], catalog['sources']['genome']['sha256'])
    accepted = [r for r in predictions if r['status'] == 'accepted']
    intact_sources = {g for g, locus in loci.items() if any(c['quality'].get('valid_orf') for c in locus['candidates'])}
    for row in accepted:
        loci[row['gene_id']]['candidates'].append(row['candidate'])
    selection = select_representatives([catalog, *donors], edges, policy=ref.DEFAULTS['policy'],
                                      min_margin=ref.DEFAULTS['min_margin'], min_support=ref.DEFAULTS['min_support'],
                                      candidate_limit=ref.DEFAULTS['candidate_limit'])
    decisions = {r['gene_id']: r for r in selection['selections'] if r['species'] == catalog['species']}
    categories, selected_categories, recovered, selected_recovered, cases = Counter(), Counter(), set(), set(), []
    selected_gff, selected_cds = [], []
    for row in accepted:
        g, c = row['gene_id'], row['candidate']
        category, exact, owners = match_category(c, *reference)
        chosen = decisions[g]['candidate_id'] == c['candidate_id']
        publishable = decisions[g]['status'] != 'ambiguous_correspondence' or c['quality'].get('rna_supported')
        categories[category] += 1
        if g in hidden and g in exact:
            recovered.add(g)
        if chosen:
            selected_categories[category] += 1
            if g in hidden and g in exact:
                selected_recovered.add(g)
            selected_gff.append(ref.candidate_gff(loci[g], c, catalog['species']))
            selected_cds.append('>' + c['source_transcript_id'] + '\n' + c['cds'] + '\n')
        cases.append(dict(gene_id=g, candidate_id=c['candidate_id'], category=category, masked=g in hidden,
                          source_intact=g in intact_sources, chosen=chosen,
                          active_full_annotation_eligible=bool(publishable), selection=decisions[g], blocks=c['blocks'],
                          alignments=row['alignments'], overlapping_reference_loci=owners))
    (output / 'selected_predictions.gff3').write_text('##gff-version 3\n' + ''.join(selected_gff))
    (output / 'selected_predictions.cds.fa').write_text(''.join(selected_cds))
    dump_compressed(output / 'predictions.json.gz', predictions)
    rescue.atomic_json(output / 'accepted_cases.json', cases)
    rescue.atomic_json(output / 'selection.json', selection)
    return dict(hidden_exon_genes=len(hidden), exact_repair_admissions=len(recovered),
                selected_exact_repairs=len(selected_recovered), accepted_categories=dict(categories),
                selected_categories=dict(selected_categories),
                selected_intact_locus_changes=sum(r['chosen'] and r['source_intact'] for r in cases),
                selected_unmasked_locus_changes=sum(r['chosen'] and not r['masked'] for r in cases),
                selected_nonexact_masked_repairs=sum(r['chosen'] and r['masked'] and r['category'] != 'exact_any_reference_isoform' for r in cases),
                full_annotation_eligible_nonexact=sum(r['active_full_annotation_eligible'] and r['category'] != 'exact_any_reference_isoform' for r in cases),
                correspondence_scope='approximate placement graph; real donor GFF/CDS/genome and actual selector; not full cohort finalizer')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('species', 'genome', 'reference-gff', 'annevo-gff', 'output'):
        parser.add_argument('--' + name, required=True, type=str if name == 'species' else Path)
    provenance = parser.add_mutually_exclusive_group(required=True)
    provenance.add_argument('--native-experiment', type=Path, help='Original registered external native experiment evidence')
    provenance.add_argument('--receipt', type=Path, help='Successful source-bound GG native GPU receipt')
    parser.add_argument('--donor-cds', action='append', required=True)
    parser.add_argument('--donor-genome', action='append', required=True)
    parser.add_argument('--donor-gff', action='append', required=True)
    parser.add_argument('--seed', type=int, action='append')
    parser.add_argument('--fraction', type=float, default=.02)
    parser.add_argument('--cpus', type=int, default=4)
    args = parser.parse_args()
    if args.seed and len(args.seed) != len(set(args.seed)):
        raise ValueError('Duplicate mask seed')
    if not 0 < args.fraction < 1 or args.cpus < 1:
        raise ValueError('Invalid mask/resource bound')
    args.output.mkdir(parents=True, exist_ok=False)
    if args.receipt:
        table = args.output / 'native_candidates.tsv'
        table.write_text('species\tgff\treceipt\n' + '\t'.join((args.species, str(args.annevo_gff.resolve()), str(args.receipt.resolve()))) + '\n')
        frozen = ann.freeze_inputs(table, {args.species: {'genome': str(args.genome), 'genetic_code': 1}})
    else:
        frozen = external_prediction(args.native_experiment, args.species, args.genome, args.annevo_gff)
    maps = {role: dict(item.split('=', 1) for item in getattr(args, 'donor_' + role)) for role in ('cds', 'genome', 'gff')}
    if len(maps['cds']) < 2 or any(set(maps[role]) != set(maps['cds']) for role in maps) or args.species in maps['cds']:
        raise ValueError('Require matching genome/GFF/CDS for at least two distinct external donor species')
    files = [args.genome, args.reference_gff, args.annevo_gff, *frozen['files'], *(Path(p) for m in maps.values() for p in m.values())]
    hashes = rescue.digest_paths(files)
    donors = [build_catalog(s, maps['cds'][s], maps['gff'][s], maps['genome'][s]) for s in sorted(maps['cds'])]
    aliases = query_locus_map(donors)
    with indexed_genome(args.genome) as genome:
        rows = bm.reference_rows(args.reference_gff, genome)
        gff, cds = bm.write_annotation(args.output / 'complete_reference', rows, genome)
        full = build_catalog(args.species, cds, gff, args.genome)
        reference = reference_index(full)
        truth = {g['gene_id']: min((c for c in g['candidates'] if c['quality']['valid_orf']), key=lambda c: c['candidate_id'])
                 for g in full['loci'] if any(c['quality']['valid_orf'] for c in g['candidates'])}
        source_ids = {g['gene_id']: g['candidates'][0]['source_gene_id'] for g in full['loci']}
        print('SOURCE_CATALOG_READY', args.species, len(full['loci']), flush=True)
        baseline, regions, proteins, placement = bm.donor_placements(args, args.output, genome)
        dump_compressed(args.output / 'baseline.checked.json.gz', baseline)
        print('PLACEMENT_READY', args.species, len(baseline), flush=True)
        native, nomination = ann.nominate(frozen, args.species, regions, proteins, genome, ref.DEFAULTS)
        checked = []
        for model in native:
            m = rescue.validate_model(model, genome, 1, ref.DEFAULTS)
            m = rescue.complete_terminals(m, genome, 1, ref.DEFAULTS, proteins[model['evidence']['donor']][model['evidence']['query']], rescue.validate_model)
            m['problems'] = sorted(set(m['problems']) | set(model['annevo_problems']))
            checked.append(rescue.check_interval(m))
        ann.confirm_transferred_structure(checked, baseline)
        dump_compressed(args.output / 'native.checked.json.gz', checked)
        rescue.atomic_json(args.output / 'nomination.json', nomination)
        print('NATIVE_NOMINATION_READY', args.species, len(checked), flush=True)
        exon_choices = {}
        for gene in full['loci']:
            g = gene['gene_id']
            if g not in truth:
                continue
            # Hide a shared internal exon, so a retained intact alternate
            # isoform cannot masquerade as a repaired masked gene.
            shared = set.intersection(*(set(tuple(b[:2]) for b in c['blocks'][1:-1]) for c in gene['candidates']))
            eligible = sorted(b for b in shared if (b[1] - b[0]) % 3)
            if eligible:
                exon_choices[g] = eligible
        masks, results = [], []
        for seed in args.seed or [20261020, 20261021]:
            hidden, exons = bm.sampled_masks(truth, exon_choices, args.fraction, seed)
            masks.append(dict(seed=seed, hidden_genes=sorted(hidden), hidden_exon_blocks={g:sorted(v) for g,v in exons.items()}))
            mgff, _ = bm.write_annotation(args.output / ('gene_' + str(seed)), rows, genome,
                                          hidden_genes={source_ids[g] for g in hidden})
            egff, ecds = bm.write_annotation(args.output / ('exon_' + str(seed)), rows, genome,
                                             hidden_blocks={source_ids[g]: v for g,v in exons.items()})
            masked = build_catalog(args.species, ecds, egff, args.genome)
            masked['annotation_spans'] = ref.annotation_ownership_spans(egff, masked)
            for condition, models in [('without_ANNEVO', baseline), ('with_ANNEVO', baseline + checked)]:
                go, eo = mgff.parent / condition, egff.parent / condition
                go.mkdir()
                eo.mkdir()
                gm = gene_assay(models, mgff, args.species, hidden, reference, go)
                em = refinement_assay(models, copy.deepcopy(masked), donors, aliases, set(exons), reference, eo)
                result = dict(seed=seed, condition=condition, gene_masking=gm, exon_masking=em)
                results.append(result)
                rescue.atomic_json(args.output / 'partial_results.json', results)
                print('CONDITION_COMPLETE', json.dumps(result), flush=True)
        control = args.output / 'unmasked_control'
        control.mkdir()
        controls = {}
        full['annotation_spans'] = ref.annotation_ownership_spans(gff, full)
        for condition, models in [('without_ANNEVO', baseline), ('with_ANNEVO', baseline + checked)]:
            d = control / condition
            d.mkdir()
            (d / 'gene').mkdir()
            (d / 'exon').mkdir()
            controls[condition] = dict(gene=gene_assay(models, gff, args.species, set(), reference, d / 'gene'),
                                      exon=refinement_assay(models, copy.deepcopy(full), donors, aliases, set(), reference, d / 'exon'))
        if rescue.digest_paths(hashes) != hashes:
            raise ValueError('Whole-genome benchmark source inputs changed')
        report = dict(schema=1, species=args.species, genome_bases=sum(genome.get_reference_length(s) for s in genome.references),
                      genome_contigs=len(genome.references), reference_genes=len(full['loci']), intact_reference_genes=len(truth),
                      all_reference_isoforms_scored=True, exon_maskable_genes=len(exon_choices),
                      input_sha256=hashes, native_identity=frozen['predictions'][args.species]['identity'], placement=placement,
                      masks=masks, results=results, unmasked_control=controls,
                      no_synthetic_RNA=True, annotation_absence_is_unknown=True,
                      scope='whole assembly admission and selector projection with real donor structure; approximate placement graph; not full cohort finalization',
                      implementation=ref.implementation(), benchmark_sha256=rescue.digest(__file__), parameters=ref.DEFAULTS)
        rescue.atomic_json(args.output / 'comparison.json', report)
        print(json.dumps({k:v for k,v in report.items() if k not in {'input_sha256','masks','implementation'}}, indent=2))


if __name__ == '__main__':
    main()

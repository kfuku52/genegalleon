#!/usr/bin/env python3
"""Repair CoGe's nested per-CDS-fragment mRNAs only with sequence evidence.

This is an explicit input-data repair, never a gene-statistics fallback. The
source files are preserved; all candidate repairs must validate before output.
"""
import argparse
import json
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq

try:
    from extend_real_test_dataset import GeneCatalog, attributes, gff_rows, read_records, sha256
except ImportError:
    from .extend_real_test_dataset import GeneCatalog, attributes, gff_rows, read_records, sha256


def repair(gff, cds, genome, species, output, audit, gene_ids=None):
    gff, cds, genome, output, audit = map(Path, (gff, cds, genome, output, audit))
    if output.exists() or audit.exists() or output.resolve() == audit.resolve():
        raise ValueError('Repair outputs must be new, distinct files')
    catalog = GeneCatalog(gff)
    sequences = read_records(cds)
    genomes = read_records(genome)
    mapping = catalog.map_cds(sequences, species)
    requested = set(Path(gene_ids).read_text().splitlines()) if gene_ids else set(mapping)
    if not requested <= mapping.keys():
        raise ValueError(f"Unknown requested CDS IDs: {sorted(requested - mapping.keys())}")
    selected_genes = {mapping[i] for i in requested}
    cds_by_gene = defaultdict(list)
    for identifier, gene in mapping.items():
        cds_by_gene[gene].append(identifier)
    rows = list(gff_rows(gff))
    groups = defaultdict(list)
    for index, row in enumerate(rows):
        roots = catalog.row_roots(row)
        if len(roots) > 1:
            raise ValueError('Shared feature across gene roots is not repairable')
        for root in roots:
            groups[root].append(index)
    replaced, omitted, changes = {}, set(), []
    for gene, indices in groups.items():
        if gene not in selected_genes:
            continue
        mrnas = {attributes(rows[i][8])['ID']: i for i in indices if rows[i][2] == 'mRNA'}
        if not any(set(attributes(rows[i][8]).get('Parent', '').split(',')) & mrnas.keys()
                   for i in mrnas.values()):
            continue
        if len(cds_by_gene[gene]) != 1:
            raise ValueError(f'Expected one CDS sequence for nested transcript gene {gene}')
        identifier = cds_by_gene[gene][0]
        feature_indices = [i for i in indices if rows[i][2] != 'gene']
        if any(rows[i][1] != 'CoGe' or rows[i][2] not in {'mRNA','CDS','exon'} for i in feature_indices):
            raise ValueError(f'Unsupported nested annotation for {gene}')
        fids = {attributes(rows[i][8]).get('coge_fid','') for i in feature_indices}
        if len(fids) != 1 or '' in fids:
            raise ValueError(f'Nested features have different CoGe IDs for {gene}')
        parents = {t: attributes(rows[i][8]).get('Parent','') for t,i in mrnas.items()}
        first = [t for t,p in parents.items() if p == gene]
        if len(first) != 1 or any(p not in set(mrnas) | {gene} for p in parents.values()):
            raise ValueError(f'Invalid nested transcript chain for {gene}')
        visited, current = set(), first[0]
        while current:
            if current in visited:
                raise ValueError(f'Cyclic transcript chain for {gene}')
            visited.add(current)
            children = [t for t,p in parents.items() if p == current]
            if len(children) > 1:
                raise ValueError(f'Branched transcripts are not repairable for {gene}')
            current = children[0] if children else None
        if visited != set(mrnas):
            raise ValueError(f'Disconnected transcript chain for {gene}')
        blocks = []
        for transcript in mrnas:
            children = [rows[i] for i in feature_indices if rows[i][2] in {'CDS','exon'}
                        and attributes(rows[i][8]).get('Parent') == transcript]
            coding = [r for r in children if r[2] == 'CDS']
            exons = [r for r in children if r[2] == 'exon']
            if len(coding) != 1 or len(exons) != 1 or coding[0][3:5] != exons[0][3:5]:
                raise ValueError(f'Expected one matching CDS/exon per fragment for {gene}')
            blocks.append(coding[0])
        if sum(rows[i][2] in {'CDS','exon'} for i in feature_indices) != 2 * len(blocks):
            raise ValueError(f'Unexpected child features for {gene}')
        if len({(r[0], r[6]) for r in blocks}) != 1 or blocks[0][6] not in {'+','-'}:
            raise ValueError(f'Inconsistent coordinates for {gene}')
        strand = blocks[0][6]
        blocks.sort(key=lambda r:int(r[3]),reverse=strand=='-')
        genomic_order = sorted((int(r[3]),int(r[4])) for r in blocks)
        if any(a < 1 or b < a for a,b in genomic_order) or any(
                genomic_order[i+1][0] <= b for i,(a,b) in enumerate(genomic_order[:-1])):
            raise ValueError(f'Overlapping or invalid CDS fragments for {gene}')
        pieces = [genomes[r[0]][1][int(r[3])-1:int(r[4])] for r in blocks]
        if strand == '-':
            pieces = [str(Seq(s).reverse_complement()) for s in pieces]
        reconstructed = ''.join(pieces).upper()
        target = sequences[identifier][1].upper()
        if len(reconstructed) != len(target) or any(a != b and b != 'N' for a,b in zip(reconstructed, target, strict=True)):
            raise ValueError(f'Genome/CDS sequence mismatch for {identifier}')
        if not any(b in 'ACGT' for b in target):
            raise ValueError(f'No resolved sequence evidence for {identifier}')
        canonical = first[0]
        row = rows[mrnas[canonical]].copy()
        row[3:5] = [str(genomic_order[0][0]),str(genomic_order[-1][1])]
        replaced[mrnas[canonical]] = row
        omitted.update(i for t,i in mrnas.items() if t != canonical)
        for i in feature_indices:
            if rows[i][2] in {'CDS','exon'}:
                row = rows[i].copy()
                fields = row[8].split(';')
                row[8] = ';'.join('Parent='+canonical if f.startswith('Parent=') else f for f in fields)
                replaced[i] = row
        changes.append(dict(gene=gene,cds_id=identifier,transcript=canonical,
                            removed_transcripts=sorted(set(mrnas)-{canonical}),cds_blocks=len(blocks),
                            cds_length=len(target),masked_bases=target.count('N'),
                            verification='equal length and all non-N CDS bases identical to genome'))
    payload = dict(inputs={kind:dict(filename=p.name,sha256=sha256(p))
                          for kind,p in [('gff',gff),('cds',cds),('genome',genome)]},repairs=changes)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as handle:
        handle.write('##gff-version 3\n')
        for i,row in enumerate(rows):
            if i not in omitted:
                handle.write('\t'.join(replaced.get(i,row))+'\n')
    payload['output_sha256'] = sha256(output)
    with audit.open('x') as handle:
        json.dump(payload,handle,indent=2)
        handle.write('\n')
    return payload


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['gff','cds','genome','species','output','audit']:
        parser.add_argument('--'+name,required=True)
    parser.add_argument('--gene-ids',help='Optional file of exact CDS IDs limiting the explicit repair scope')
    result=repair(**vars(parser.parse_args()))
    print(f"Repaired {len(result['repairs'])} sequence-validated CoGe transcript chains")


if __name__ == '__main__':
    main()

"""Regression contract for the checked-in biological test input."""

import hashlib
import json
from pathlib import Path

from workflow.support.extend_real_test_dataset import GeneCatalog, gff_rows, read_records, sha256

INPUT = Path(__file__).resolve().parents[2] / "workspace/input"


def test_real_aha_dataset_preserves_seed_genes_and_unrelated_inputs():
    manifest = json.loads((INPUT / "dataset_manifest/real_neighborhoods.json").read_text())
    assert manifest["neighbors"] == 20
    for relative, digest in manifest["preserved_files"].items():
        assert sha256(INPUT / relative) == digest, relative
    for species in manifest["species"]:
        files = {kind: INPUT / kind / entry["filename"]
                 for kind, entry in species["output_files"].items()}
        cds = read_records(files["species_cds"])
        assert set(species["retained_ids"]) <= cds.keys()
        assert not set(species["removed_dummy_ids"]) & cds.keys()
        for identifier, digest in species["retained_cds_sha256"].items():
            assert hashlib.sha256(cds[identifier][1].encode()).hexdigest() == digest
        for kind, path in files.items():
            assert sha256(path) == species["output_files"][kind]["sha256"]


def test_real_aha_neighborhoods_are_contiguous_and_annotations_fit():
    manifest = json.loads((INPUT / "dataset_manifest/real_neighborhoods.json").read_text())
    anchors = set((INPUT / "dataset_manifest/aha_anchors.txt").read_text().splitlines())
    assert {entry["anchor"] for entry in manifest["coverage"]} == anchors
    for species in manifest["species"]:
        name = species["species"]
        files = {kind: INPUT / kind / entry["filename"]
                 for kind, entry in species["output_files"].items()}
        genome = read_records(files["species_genome"])
        catalog = GeneCatalog(files["species_gff"])
        mapping = catalog.map_cds(read_records(files["species_cds"]), name)
        for row in gff_rows(files["species_gff"]):
            assert "gg_dummy_" not in row[8]
            assert 1 <= int(row[3]) <= int(row[4]) <= len(genome[row[0]][1])
        for entry in manifest["coverage"]:
            if entry["species"] != name:
                continue
            gene = mapping[entry["anchor"]]
            chromosome = catalog.genes[gene][0]
            for neighbor in entry["left_genes"] + entry["right_genes"]:
                assert catalog.genes[neighbor][0] == chromosome
            assert entry["left_count"] == len(entry["left_genes"]) <= 20
            assert entry["right_count"] == len(entry["right_genes"]) <= 20
        for window in manifest["windows"]:
            if window["species"] == name:
                seq = genome[window["window_id"]][1]
                assert len(seq) == window["end"] - window["start"] + 1
                assert hashlib.sha256(seq.encode()).hexdigest() == window["sequence_sha256"]


def test_all_aha_structures_reconstruct_the_input_cds():
    """Catch wrong isoforms/fragmented models even when IDs resolve correctly."""
    from Bio.Seq import Seq

    from workflow.support.gff2genestat import process_single_gff
    manifest=json.loads((INPUT/'dataset_manifest/real_neighborhoods.json').read_text())
    anchors=set((INPUT/'dataset_manifest/aha_anchors.txt').read_text().splitlines())
    checked=set()
    columns=['gene_id','feature_size','num_intron','feature_blocks','strand','chromosome']
    gff_columns=['sequence','source','feature','start','end','score','strand','phase','attributes']
    for species in manifest['species']:
        files={kind:INPUT/kind/entry['filename'] for kind,entry in species['output_files'].items()}
        cds=read_records(files['species_cds'])
        genome=read_records(files['species_genome'])
        selected=sorted(anchors & cds.keys())
        traits=process_single_gff(files['species_gff'].name,str(files['species_gff'].parent),
                                  selected,'CDS','longest',gff_columns,columns)
        for row in traits.itertuples(index=False):
            blocks=[tuple(map(int,b.split('-'))) for b in row.feature_blocks.split(';')]
            pieces=[genome[row.chromosome][1][a-1:b] for a,b in blocks]
            if row.strand=='-':
                pieces=[str(Seq(s).reverse_complement()) for s in pieces]
            sequence=''.join(pieces).upper()
            target=cds[row.gene_id][1].upper()
            assert len(sequence)==len(target)==row.feature_size, row.gene_id
            assert all(a==b or b=='N' for a,b in zip(sequence, target, strict=True)), row.gene_id
            if row.gene_id.startswith('Cephalotus_'):
                assert row.num_intron>0
            checked.add(row.gene_id)
    assert checked==anchors and len(checked)==60

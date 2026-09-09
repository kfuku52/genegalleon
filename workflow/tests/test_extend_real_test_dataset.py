import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from workflow.support.extend_real_test_dataset import GeneCatalog, build_species, read_records


def make_inputs(tmp_path):
    source, seed, out = [tmp_path / name for name in ["source", "seed", "out"]]
    for root in [source, seed]:
        for kind in ["species_cds", "species_gff", "species_genome"]:
            (root / kind).mkdir(parents=True)
    rows, cds = [], []
    # g0/g44 are scaffold-edge cases; the second contig must never be bridged.
    for i in range(50):
        chrom = "chr1" if i < 45 else "chr2"
        start, end = i * 100 + 11, i * 100 + 40
        strand = "-" if i == 22 else "+"
        for feature, attr in [("gene", f"ID=g{i}"),
                              ("mRNA", f"ID=t{i};Parent=g{i}"),
                              ("CDS", f"ID=c{i};Parent=t{i}")]:
            rows.append(f"{chrom}\treal\t{feature}\t{start}\t{end}\t.\t{strand}\t0\t{attr}\n")
        cds.append(f">Species_a_t{i}\n" + "ACG" * 10 + "\n")
    (source / "species_gff/Species_a.gff").write_text("".join(rows))
    (source / "species_cds/Species_a.fa").write_text("".join(cds))
    (source / "species_genome/Species_a.fa").write_text(
        ">chr1\n" + "ACGT" * 1250 + "\n>chr2\n" + "TGCA" * 1250 + "\n")
    # Existing non-AHA gene on another contig plus one marked artificial gene.
    (seed / "species_cds/Species_a.fa").write_text(
        cds[22] + cds[48] + ">Species_a_fake\nATGAAA\n")
    (seed / "species_gff/Species_a.gff").write_text(
        "chr1\tCoGe\tCDS\t1\t6\t.\t+\t.\tAlias=fake;coge_fid=gg_dummy_fake\n")
    (seed / "species_genome/Species_a.fa").write_text(">old\nACGT\n")
    (seed / "unrelated.txt").write_text("another family's input\n")
    shutil.copytree(seed, out)
    return source, seed, out


def test_complete_neighborhoods_preserve_other_genes_and_source_coordinates(tmp_path):
    source, seed, out = make_inputs(tmp_path)
    seed_bytes = (seed / "species_cds/Species_a.fa").read_bytes()
    summary, coverage, windows = build_species(
        source, seed, out, "Species_a", {"Species_a_t22"}, 20, 5)
    records = read_records(out / "species_cds/Species_a.fa")
    assert set(records) == {f"Species_a_t{i}" for i in range(2, 43)} | {"Species_a_t48"}
    assert summary["retained_ids"] == ["Species_a_t22", "Species_a_t48"]
    assert coverage[0]["left_count"] == coverage[0]["right_count"] == 20
    assert coverage[0]["strand"] == "-"
    assert len([w for w in windows if w["source_seqid"] == "chr1"]) == 1
    assert (out / "unrelated.txt").read_bytes() == (seed / "unrelated.txt").read_bytes()
    assert (seed / "species_cds/Species_a.fa").read_bytes() == seed_bytes
    assert "Species_a_fake" in read_records(out / "dataset_manifest/synthetic_synteny/Species_a.fa")
    originals = set((source / "species_gff/Species_a.gff").read_text().splitlines())
    converted = []
    for line in (out / "species_gff/Species_a.gff").read_text().splitlines()[1:]:
        row = line.split("\t")
        chrom, coords = row[0].rsplit(":", 1)
        start, end = map(int, coords.split("-"))
        assert 1 <= int(row[3]) <= int(row[4]) <= end - start + 1
        row[0] = chrom
        row[3:5] = [str(int(row[3]) + start - 1), str(int(row[4]) + start - 1)]
        converted.append("\t".join(row))
    assert set(converted) <= originals
    assert len(converted) == 42 * 3
    genome = read_records(out / "species_genome/Species_a.fa")
    full = read_records(source / "species_genome/Species_a.fa")
    for window in windows:
        assert genome[window["window_id"]][1] == full[window["source_seqid"]][1][window["start"] - 1:window["end"]]


def test_scaffold_edges_are_reported_without_fabricating_neighbors(tmp_path):
    source, seed, out = make_inputs(tmp_path)
    with (seed / "species_cds/Species_a.fa").open("a") as handle:
        handle.write(">Species_a_t0\n" + "ACG" * 10 + "\n")
    _, coverage, _ = build_species(source, seed, out, "Species_a", {"Species_a_t0"}, 20, 5)
    assert (coverage[0]["left_count"], coverage[0]["right_count"]) == (0, 20)
    assert coverage[0]["right_genes"] == [f"g{i}" for i in range(1, 21)]


@pytest.mark.parametrize("defect", ["missing", "changed"])
def test_existing_non_anchor_gene_cannot_silently_drop_or_change(tmp_path, defect):
    source, seed, out = make_inputs(tmp_path)
    path = source / "species_cds/Species_a.fa"
    text = path.read_text()
    old = ">Species_a_t48\n" + "ACG" * 10 + "\n"
    path.write_text(text.replace(old, "" if defect == "missing" else old.replace("ACG", "AAA")))
    with pytest.raises(ValueError, match="Existing"):
        build_species(source, seed, out, "Species_a", {"Species_a_t22"}, 20, 5)


def test_boundary_overlaps_keep_whole_gene_models(tmp_path):
    source, seed, out = make_inputs(tmp_path)
    # Padding reaches the preceding gene, so its full model must be retained.
    build_species(source, seed, out, "Species_a", {"Species_a_t22"}, 1, 80)
    assert "Species_a_t20" in read_records(out / "species_cds/Species_a.fa")
    catalog = GeneCatalog(out / "species_gff/Species_a.gff")
    assert catalog.genes["g20"][2] - catalog.genes["g20"][1] + 1 == 30


def test_neighbor_order_uses_cds_bounds_when_utrs_overlap(tmp_path):
    source, seed, out = make_inputs(tmp_path)
    gff = source / "species_gff/Species_a.gff"
    gff.write_text(gff.read_text().replace(
        "gene\t2111\t2140\t.\t+\t0\tID=g21",
        "gene\t1991\t2140\t.\t+\t0\tID=g21"))
    _, coverage, _ = build_species(source, seed, out, "Species_a", {"Species_a_t22"}, 1, 5)
    assert coverage[0]["left_genes"] == ["g21"]
    assert coverage[0]["right_genes"] == ["g23"]


def test_minimal_builder_does_not_clip_boundary_gene_models(tmp_path):
    from workflow.support.build_minimal_test_dataset import EffectiveWindow, write_shifted_gff

    source, _, _ = make_inputs(tmp_path)
    window = EffectiveWindow(seqid="chr1", start=2115, end=2245,
                             window_id="chr1:2115-2245", cores={"t22"})
    output = tmp_path / "minimal.gff"
    write_shifted_gff(source / "species_gff/Species_a.gff", {"chr1": [window]}, output)
    catalog = GeneCatalog(output)
    assert set(catalog.genes) == {"g22"}
    assert catalog.genes["g22"][1:3] == (97, 126)
    assert len(output.read_text().splitlines()) == 4  # Header and three complete features.


def test_cli_refuses_existing_output_and_preserves_unrelated_inputs(tmp_path):
    source, seed, out = make_inputs(tmp_path)
    anchors = tmp_path / "anchors.txt"
    anchors.write_text("Species_a_t22\n")
    script = Path(__file__).resolve().parents[1] / "support/extend_real_test_dataset.py"
    cmd = [sys.executable, str(script), "--source-pg", str(source), "--seed-pg", str(seed),
           "--out-pg", str(out), "--anchor-ids", str(anchors)]
    failed = subprocess.run(cmd, capture_output=True, text=True)
    assert failed.returncode != 0
    assert "new directory" in failed.stderr
    shutil.rmtree(out)
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    manifest = json.loads((out / "dataset_manifest/real_neighborhoods.json").read_text())
    assert manifest["neighbors"] == 20
    assert (out / "unrelated.txt").read_bytes() == (seed / "unrelated.txt").read_bytes()


def test_adds_missing_gff_and_genome_without_replacing_seed_cds(tmp_path):
    source,seed,out=make_inputs(tmp_path)
    for kind in ['species_gff','species_genome']:
        (seed/kind/'Species_a.gff' if kind=='species_gff' else seed/kind/'Species_a.fa').unlink()
    cds=seed/'species_cds/Species_a.fa'
    cds.write_text(cds.read_text().split('>Species_a_fake')[0])
    summary,_,_=build_species(source,seed,out,'Species_a',{'Species_a_t22'},1,5)
    assert set(summary['seed_files'])=={'species_cds'}
    assert (out/'species_gff/Species_a.gff').exists()
    assert (out/'species_genome/Species_a.fa').exists()
    assert read_records(cds).items() <= read_records(out/'species_cds/Species_a.fa').items()

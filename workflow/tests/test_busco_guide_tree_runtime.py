"""Real compiled k-mer and RapidNJ tests on bounded synthetic inputs."""
import csv
import json
import math
import random
import shutil
import struct
import subprocess

import pytest
from Bio import Phylo

from workflow.support import busco_guide_tree as guide
from workflow.tests.test_busco_guide_tree import make_archive


def run_kernel(args):
    result = subprocess.run(["gg-kmer-distance", *map(str, args)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def reference_hashes(sequence, k, size):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    mask = (1 << 64) - 1
    values = set()
    for start in range(len(sequence)-k+1):
        fragment = sequence[start:start+k]
        if not all(c in alphabet for c in fragment):
            continue
        x = 0
        for c in fragment:
            x = x*20 + alphabet.index(c)
        x = (x + 0x9e3779b97f4a7c15) & mask
        x = ((x ^ (x >> 30)) * 0xbf58476d1ce4e5b9) & mask
        x = ((x ^ (x >> 27)) * 0x94d049bb133111eb) & mask
        values.add(x ^ (x >> 31))
    return sorted(values)[:size]


def decode(path):
    raw = path.read_bytes()
    assert raw[:8] == b"GGKMER01"
    k, size, markers = struct.unpack_from("<QQQ", raw, 8)
    position, data = 32, []
    for _ in range(markers):
        count, = struct.unpack_from("<Q", raw, position)
        position += 8
        data.append(list(struct.unpack_from("<" + "Q"*count, raw, position)))
        position += 8*count
    assert position == len(raw)
    return k, size, data


@pytest.mark.parametrize("size", [16, 256])
def test_kernel_matches_independent_hash_and_union_distance(tmp_path, size):
    sequences = ["ACDEFGHIKLMNPQRSTVWY"*4, "ACDEFGHIKLXXPQRSTVWY"*4, "YYYYYYYYYY"]
    paths = []
    for i, seq in enumerate(sequences):
        source, binary = tmp_path / f"{i}.txt", tmp_path / f"{i}.bin"
        source.write_text("2\n" + seq + "\n" + (seq if i != 2 else "") + "\n")
        run_kernel(["sketch", source, binary, 5, size])
        assert decode(binary)[2][0] == reference_hashes(seq, 5, size)
        paths.append(binary)
    listing = tmp_path / "list.txt"
    listing.write_text("\n".join(map(str, paths)) + "\n")
    outputs = []
    for cpus in (1, 3):
        output = tmp_path / f"distances{cpus}.tsv"
        run_kernel(["compare", listing, output, cpus])
        outputs.append(output.read_bytes())
    assert outputs[0] == outputs[1]
    with (tmp_path / "distances1.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        a, b = [set(reference_hashes(sequences[int(row[key])], 5, size)) for key in ("i", "j")]
        sample = set(sorted(a | b)[:size])
        similarity = len(a & b & sample)/len(sample)
        expected = min(1, -math.log(2*similarity/(1+similarity))/5) if similarity else 1
        assert float(row["distance"]) == pytest.approx(expected)
    # Missing markers do not count; saturated present proteins still do.
    assert int(rows[1]["shared"]) == 1 and int(rows[1]["saturated"]) == 1


def test_kernel_rejects_corrupt_binary_and_invalid_arguments(tmp_path):
    bad = tmp_path / "bad.bin"
    bad.write_bytes(b"GGKMER01")
    listing = tmp_path / "list.txt"
    listing.write_text(str(bad) + "\n")
    result = subprocess.run(["gg-kmer-distance", "compare", str(listing), str(tmp_path / "out"), "1"], capture_output=True, text=True)
    assert result.returncode and "Truncated" in result.stderr


def fixture(tmp_path):
    rng = random.Random(417)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    prototypes = {f"BUSCO{g:04d}": "M" + "".join(rng.choices(alphabet, k=249)) for g in range(40)}
    names = [f"Plant_species{i:03d}" for i in range(6)]
    for i, name in enumerate(names):
        sequences = {}
        for marker, sequence in prototypes.items():
            common = random.Random(1000*int(marker[5:]) + i//2)
            dna = list(sequence)
            for index in common.sample(range(1,250), 30):
                dna[index] = common.choice(alphabet)
            private = random.Random(10000*int(marker[5:]) + i)
            for index in private.sample(range(1,250), 2):
                dna[index] = private.choice(alphabet)
            sequences[marker] = "".join(dna)
        make_archive(tmp_path, name, sequences)
    return names


def build_args(tmp_path, output="guide"):
    return guide.parser().parse_args(["build", "--cds-dir", str(tmp_path / "cds"), "--full-dir", str(tmp_path / "full"),
                                      "--short-dir", str(tmp_path / "short"), "--output", str(tmp_path / output),
                                      "--cache", str(tmp_path / "cache"), "--markers", "40", "--minimum-shared", "10",
                                      "--nearest", "1", "--cpus", "2"])


def test_real_rapidnj_guide_neighbors_cache_and_freezing(tmp_path):
    assert shutil.which("rapidnj") and shutil.which("gg-kmer-distance")
    names = fixture(tmp_path)
    args = build_args(tmp_path)
    first = guide.build(args)
    stability = json.loads((args.output / "stability.json").read_text())
    for i, name in enumerate(names):
        assert stability[name]["nearest"] == [names[i ^ 1]]
        assert stability[name]["stable"]
    tree = Phylo.read(args.output / "guide_tree.nwk", "newick")
    assert {tip.name for tip in tree.get_terminals()} == set(names)
    assert all((c.branch_length or 0) >= 0 for c in tree.find_clades())
    assert first["performance"]["sketches_reused"] == 0
    second = guide.build(build_args(tmp_path, "guide2"))
    assert second["performance"]["sketches_reused"] == len(names)
    assert second["performance"]["distances_reused"]
    assert (tmp_path / "guide/guide_tree.nwk").read_bytes() == (tmp_path / "guide2/guide_tree.nwk").read_bytes()
    args.k = 6
    with pytest.raises(ValueError, match="Frozen BUSCO guide differs"):
        guide.build(args)


def test_changed_sources_and_too_few_shared_fail_without_receipt(tmp_path):
    fixture(tmp_path)
    args = build_args(tmp_path)
    args.minimum_shared = 100
    with pytest.raises(ValueError, match="Too few high-occupancy"):
        guide.build(args)
    assert not (args.output / "receipt.json").exists()
    assert list(args.output.glob(".failed*"))


def test_changed_archive_is_detected_during_execution(tmp_path, monkeypatch):
    names = fixture(tmp_path)
    original = guide.execute
    def replace_after_compare(command, root, label):
        result = original(command, root, label)
        if label == "compare":
            archive = tmp_path / "full/single_copy" / (names[0] + ".json.gz")
            archive.write_bytes(archive.read_bytes() + b"changed")
        return result
    monkeypatch.setattr(guide, "execute", replace_after_compare)
    with pytest.raises(OSError, match="File changed"):
        guide.build(build_args(tmp_path))
    assert not (tmp_path / "guide/receipt.json").exists()


def test_zero_information_tree_fails_without_taxonomy_fallback(tmp_path):
    for i in range(3):
        make_archive(tmp_path, f"Plant_species{i:03d}", {f"BUSCO{g:04d}": "MDEKAAA" for g in range(40)})
    args = build_args(tmp_path)
    with pytest.raises(ValueError, match="No informative guide-tree"):
        guide.build(args)
    assert not (args.output / "receipt.json").exists()


def test_corrupt_sketch_and_matrix_are_recomputed(tmp_path):
    names = fixture(tmp_path)
    args = build_args(tmp_path)
    guide.build(args)
    matrix_metadata = next(args.cache.glob("*.matrix.json"))
    key = json.loads(matrix_metadata.read_text())["key"]
    (args.cache / (key+".tsv")).write_text("corrupt matrix\n")
    binary = next(args.cache.glob("*.bin"))
    binary.write_bytes(b"corrupt sketch")
    second = guide.build(build_args(tmp_path, "recomputed"))
    assert second["performance"]["sketches_reused"] == len(names)-1
    assert not second["performance"]["distances_reused"]
    assert (tmp_path / "guide/guide_tree.nwk").read_bytes() == (tmp_path / "recomputed/guide_tree.nwk").read_bytes()


def test_all_saturated_markers_fail_instead_of_arbitrary_neighbors(tmp_path):
    for i, sequence in enumerate(("MAAAAAA", "CCCCCCC", "DDDDDDD")):
        make_archive(tmp_path, f"Plant_species{i:03d}", {f"BUSCO{g:04d}": sequence for g in range(40)})
    args = build_args(tmp_path)
    with pytest.raises(ValueError, match="All usable BUSCO marker pairs are k-mer saturated"):
        guide.build(args)
    assert not (args.output / "receipt.json").exists()


def test_guide_connects_real_reference_planning_and_bounded_alternatives(tmp_path):
    from workflow.support import rescue_gene_models as rescue
    names = fixture(tmp_path)
    args = build_args(tmp_path)
    guide.build(args)
    # Controlled panel disagreement: its extra candidate must survive planning.
    stability_path = args.output / "stability.json"
    stability = json.loads(stability_path.read_text())
    stability[names[0]]["panel_0"] = [names[2]]
    stability[names[0]]["stable"] = False
    stability_path.write_text(json.dumps(stability))
    receipt_path = args.output / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["files"]["stability.json"] = guide.digest(stability_path)
    receipt_path.write_text(json.dumps(receipt))
    for folder in ("gff", "genome"):
        (tmp_path / folder).mkdir()
    for name in names:
        (tmp_path / "genome" / (name+".fa")).write_text(">chr1\nATGAAATAA\n")
        (tmp_path / "gff" / (name+".gff3")).write_text(
            f"chr1\tfixture\tCDS\t1\t9\t.\t+\t0\tID={name}_gene\n")
    options = rescue.parser().parse_args(["plan", "--cds-dir", str(tmp_path / "cds"),
        "--gff-dir", str(tmp_path / "gff"), "--genome-dir", str(tmp_path / "genome"),
        "--busco-dir", str(tmp_path / "short"), "--tree", str(args.output / "guide_tree.nwk"),
        "--guide-tree-receipt", str(receipt_path), "--output", str(tmp_path / "rescue"),
        "--nearest-references", "1"])
    plan = rescue.build_plan(options)
    assert plan["nearest_references"][names[0]] == [names[1], names[2]]
    assert plan["nearest_references"][names[1]] == [names[0]]
    assert plan["request"]["tree_metric"] == "patristic_distance"
    assert len(plan["common_references"]) == 5
    assert plan["guide_tree_stability"][names[0]]["stable"] is False
    # A different CDS set, even with the same species names, is inadmissible.
    cds = tmp_path / "cds" / (names[0]+".fa")
    cds.write_text(cds.read_text()+"\n")
    options.output = tmp_path / "wrong-input-rescue"
    with pytest.raises(ValueError, match="Frozen BUSCO guide inputs"):
        rescue.build_plan(options)

import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "support" / "shorten_fasta_newick_names.sh"


def _run_script(
    in_fasta: Path,
    out_fasta: Path,
    in_newick: Path,
    out_newick: Path,
    max_length: int,
) -> None:
    subprocess.run(
        [
            "bash",
            str(SCRIPT),
            str(in_fasta),
            str(out_fasta),
            str(in_newick),
            str(out_newick),
            str(max_length),
        ],
        check=True,
        cwd=str(REPO_ROOT),
        text=True,
        capture_output=True,
    )


def test_shorten_fasta_newick_names_rewrites_long_ids(tmp_path: Path):
    in_fasta = tmp_path / "in.fa"
    out_fasta = tmp_path / "out.fa"
    in_newick = tmp_path / "in.nwk"
    out_newick = tmp_path / "out.nwk"

    in_fasta.write_text(
        (
            ">long.name+alpha123\n"
            "ATGC\n"
            ">short\n"
            "ATGC\n"
            ">long.name+beta456\n"
            "ATGC\n"
        ),
        encoding="utf-8",
    )
    in_newick.write_text(
        "(long.name+alpha123:0.1,short:0.2,long.name+beta456:0.3);\n",
        encoding="utf-8",
    )

    _run_script(in_fasta, out_fasta, in_newick, out_newick, max_length=10)

    fasta_out = out_fasta.read_text(encoding="utf-8")
    newick_out = out_newick.read_text(encoding="utf-8")

    # NR-based suffixing: long.name+alpha123 -> long.nam_1, long.name+beta456 -> long.nam_2
    assert ">long.nam_1" in fasta_out
    assert ">long.nam_2" in fasta_out
    assert "long.name+alpha123" not in fasta_out
    assert "long.name+beta456" not in fasta_out
    assert "long.nam_1:0.1" in newick_out
    assert "long.nam_2:0.3" in newick_out
    assert "long.name+alpha123" not in newick_out
    assert "long.name+beta456" not in newick_out


def test_shorten_fasta_newick_names_keeps_inputs_when_no_long_ids(tmp_path: Path):
    in_fasta = tmp_path / "in.fa"
    out_fasta = tmp_path / "out.fa"
    in_newick = tmp_path / "in.nwk"
    out_newick = tmp_path / "out.nwk"

    fasta_text = ">sp1\nATGC\n>sp2\nATGC\n"
    newick_text = "(sp1:0.1,sp2:0.2);\n"
    in_fasta.write_text(fasta_text, encoding="utf-8")
    in_newick.write_text(newick_text, encoding="utf-8")

    _run_script(in_fasta, out_fasta, in_newick, out_newick, max_length=10)

    assert out_fasta.read_text(encoding="utf-8") == fasta_text
    assert out_newick.read_text(encoding="utf-8") == newick_text


def test_shortening_quoted_newick_tips_uses_the_same_mapping(tmp_path):
    from ete4 import Tree
    fasta, tree = tmp_path / "in.fa", tmp_path / "in.nwk"
    output, output_tree = tmp_path / "out.fa", tmp_path / "out.nwk"
    fasta.write_text(">abcdefghijk\nATG\n>short\nATG\n")
    tree.write_text("('abcdefghijk':1,short:2);\n")
    _run_script(fasta, output, tree, output_tree, 8)
    tips = set(Tree(output_tree.read_text(), parser=1).leaf_names())
    assert tips == {line[1:] for line in output.read_text().splitlines() if line.startswith(">")}
    assert tips == {"abcdef_1", "short"}


@pytest.mark.parametrize('maximum', [8, 20])
def test_shortening_preserves_crlf_sequence_bytes_and_missing_final_newline(tmp_path, maximum):
    source, tree = tmp_path / 'in.fa', tmp_path / 'in.nwk'
    output, output_tree = tmp_path / 'out.fa', tmp_path / 'out.nwk'
    data = b'>abcdefghijk\r\nAc-G\r\n>short\r\nA?TG'
    source.write_bytes(data)
    tree.write_text("('abcdefghijk':1,short:2);")
    _run_script(source, output, tree, output_tree, maximum)
    expected = data.replace(b'abcdefghijk', b'abcdef_1') if maximum == 8 else data
    assert output.read_bytes() == expected
    assert source.read_bytes() == data


def test_shortening_failure_preserves_the_existing_fasta_tree_pair(tmp_path, monkeypatch):
    import importlib.util
    path = SCRIPT.with_suffix('.py')
    spec = importlib.util.spec_from_file_location('shorten_adapter', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source, tree = tmp_path / 'in.fa', tmp_path / 'in.nwk'
    output, output_tree = tmp_path / 'out.fa', tmp_path / 'out.nwk'
    source.write_text('>abcdefghijk\nATG\n>short\nATG\n')
    tree.write_text('(abcdefghijk:1,short:2);')
    output.write_bytes(b'previous fasta')
    output_tree.write_bytes(b'previous tree')
    def fail(_args):
        raise ValueError('tree serialization failure')
    monkeypatch.setattr(module, 'rename_main', fail)
    with pytest.raises(ValueError, match='tree serialization failure'):
        module.shorten(source, output, tree, output_tree, 8)
    assert output.read_bytes() == b'previous fasta'
    assert output_tree.read_bytes() == b'previous tree'


def test_shorten_fasta_newick_names_avoids_collisions_with_existing_ids(tmp_path: Path):
    in_fasta = tmp_path / "in.fa"
    out_fasta = tmp_path / "out.fa"
    in_newick = tmp_path / "in.nwk"
    out_newick = tmp_path / "out.nwk"

    in_fasta.write_text(
        (
            ">long.name+alpha123\n"
            "ATGC\n"
            ">long.nam_1\n"
            "ATGC\n"
            ">long.name+beta456\n"
            "ATGC\n"
        ),
        encoding="utf-8",
    )
    in_newick.write_text(
        "(long.name+alpha123:0.1,long.nam_1:0.2,long.name+beta456:0.3);\n",
        encoding="utf-8",
    )

    _run_script(in_fasta, out_fasta, in_newick, out_newick, max_length=10)

    fasta_out = out_fasta.read_text(encoding="utf-8")
    newick_out = out_newick.read_text(encoding="utf-8")

    # "long.nam_1" already exists, so generated IDs must skip to non-colliding values.
    assert ">long.nam_1" in fasta_out
    assert ">long.nam_2" in fasta_out
    assert ">long.nam_3" in fasta_out
    assert "long.name+alpha123" not in fasta_out
    assert "long.name+beta456" not in fasta_out
    assert "(long.nam_2:0.1,long.nam_1:0.2,long.nam_3:0.3);" in newick_out

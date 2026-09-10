"""Native CAFE adapter contracts; evolutionary likelihoods are not reimplemented."""
import json
from pathlib import Path
from types import SimpleNamespace
import pytest

from workflow.support import cafe_branch_specificity as m

ASR = """#nexus
BEGIN TREES;
 TREE F1 = ((A<1>*_9:1,B<2>_3:1)<5>_4:1,C<3>_4:2)<4>_4;
END;
"""


@pytest.fixture
def tree_file(tmp_path):
    path = tmp_path / "Gamma_asr.tre"
    path.write_text(ASR)
    return path


def fit(score, lambdas=(.1,), directory="native"):
    return m.Fit(score, list(lambdas), 10, directory)


def test_native_branch_mapping_uses_ids_not_substrings_or_species_order(tree_file):
    t = m.read_cafe_tree(tree_file, "<5>")
    assert t.newick == "((A:1,B:1):1,C:2);\n"
    assert t.lambda_newick == "((A:1,B:1):2,C:1):1;\n"
    assert t.root_label == "<4>"
    assert t.parent_label == "<4>"
    assert t.tip_labels == {"A": "A<1>", "B": "B<2>", "C": "C<3>"}
    assert m.read_cafe_tree(tree_file, "A<1>").lambda_newick == "((A:2,B:1):1,C:1):1;\n"
    for invalid in ("A", "<1>", "<4>", "1"):
        with pytest.raises(ValueError, match="Target"):
            m.read_cafe_tree(tree_file, invalid)


def test_original_durations_are_matched_by_clade_not_order(tree_file, tmp_path):
    dated = tmp_path / "dated.nwk"
    dated.write_text("(C:2,(B:1,A:1.000001):1);")
    tree = m.read_cafe_tree(tree_file, "<5>", dated)
    assert "A:1.0000009999999999" in tree.newick
    dated.write_text("((A:1,C:1):1,B:2);")
    with pytest.raises(ValueError, match="topology"):
        m.read_cafe_tree(tree_file, "<5>", dated)
    dated.write_text("(C:2,(B:1,A:10):1);")
    with pytest.raises(ValueError, match="durations"):
        m.read_cafe_tree(tree_file, "<5>", dated)


@pytest.mark.parametrize("replace", [
    ("A<1>*_9:1", "A<1>*_9:0"),
    ("B<2>_3", "A<1>_3"),
    ("C<3>_4", "C_4"),
])
def test_invalid_asr_is_rejected(tree_file, replace):
    tree_file.write_text(ASR.replace(*replace))
    with pytest.raises(ValueError):
        m.read_cafe_tree(tree_file, "<5>")


def test_lr_uses_native_nll_sign_and_rejects_unresolved_nested_fit():
    assert m.lr_statistic(fit(12), fit(9, (.1, .8))) == 6
    assert m.lr_statistic(fit(12), fit(12 + 1e-7, (.1, .1))) == 0
    with pytest.raises(m.NativeCafeError, match="worse"):
        m.lr_statistic(fit(12), fit(13, (.1, .1)))


def test_bootstrap_resolution_depends_on_replicates_and_includes_ties():
    assert m.bootstrap_pvalue(10, [0] * 99) == .01
    assert m.bootstrap_pvalue(10, [0, 10, 11]) == .75
    assert m.bootstrap_pvalue(0, [0, 1]) == 1
    for invalid in ([], [float("nan")], [-1]):
        with pytest.raises(ValueError):
            m.bootstrap_pvalue(1, invalid)


def make_fit_outputs(directory, score="10", lambdas="0.1, 0.2", iterations=20):
    (directory / "results").mkdir(parents=True, exist_ok=True)
    (directory / "cafe.log").write_text(f"Completed {iterations} iterations\nFinal -lnL: {score}\n")
    (directory / "results" / "Base_lambda_per_family.txt").write_text(f"family\t{lambdas}\n")


def test_native_fit_parser_rejects_nonfinite_truncated_or_capped_results(tmp_path):
    make_fit_outputs(tmp_path)
    assert m.parse_fit(tmp_path, 2, 100).lambdas == [.1, .2]
    for kwargs in ({"score": "nan"}, {"lambdas": "0.1"}, {"iterations": 100}, {"lambdas": "0.1, -0.2"}):
        make_fit_outputs(tmp_path, **kwargs)
        with pytest.raises(m.NativeCafeError):
            m.parse_fit(tmp_path, 2, 100)
    (tmp_path / "cafe.log").write_text("interrupted")
    with pytest.raises(m.NativeCafeError, match="Incomplete"):
        m.parse_fit(tmp_path, 2, 100)


def test_both_models_are_refit_for_every_bootstrap_dataset(tmp_path):
    class Engine:
        def __init__(self):
            self.calls = []
        def fit(self, counts, alternative):
            self.calls.append((tuple(counts), alternative))
            nll = {1: (12, 9), 2: (11, 10), 3: (13, 9)}[counts[0]][alternative]
            return fit(nll, (.1, .2) if alternative else (.1,))
        def null_root(self, counts, null):
            return 4, "null_reconstruction"
        def simulate(self, counts, null, root, replicates):
            assert root == 4 and replicates == 2
            return [[2], [3]], "simulation"
    engine = Engine()
    result = m.compare_family(engine, [1], 2, tmp_path)
    assert engine.calls == [((i,), alt) for i in (1, 2, 3) for alt in (False, True)]
    assert result["p_value"] == pytest.approx(2 / 3)
    assert result["bootstrap_completed"] == 2
    assert (tmp_path / "bootstrap.tsv").read_text().count("\n") == 3


def test_zero_statistic_has_exact_p_one_without_fake_simulations(tmp_path):
    class Engine:
        def fit(self, counts, alternative):
            return fit(10, (.1, .1) if alternative else (.1,))
    result = m.compare_family(Engine(), [1], 999, tmp_path)
    assert result["p_value"] == 1
    assert result["p_method"] == "zero_statistic"
    assert result["bootstrap_completed"] == 0


@pytest.fixture
def native_stub(tmp_path):
    path = tmp_path / "cafe5"
    path.write_text("""#!/usr/bin/env python3
import pathlib, sys
args = sys.argv[1:]
out = pathlib.Path(args[args.index('-o') + 1]); out.mkdir(parents=True)
assert '-b' in args and '-i' in args and '-t' in args and '-c' in args
alt = '-y' in args
(out / 'Base_lambda_per_family.txt').write_text('family\\t' + ('0.1, 0.2' if alt else '0.1') + '\\n')
print('Completed 20 iterations')
print('Final -lnL: ' + ('9' if alt else '10'))
""")
    path.chmod(0o755)
    return path


def test_native_commands_restarts_and_content_verified_resume(tree_file, native_stub, tmp_path):
    tree = m.read_cafe_tree(tree_file, "<5>")
    engine = m.NativeCafe(tree, tmp_path / "work", str(native_stub), 2, 100, 1, 10)
    first = engine.fit([9, 3, 4], True)
    commands = list(engine.cache.glob("*/attempt_*/command.json"))
    assert len(commands) == 2
    for path in commands:
        argv = json.loads(path.read_text())["argv"]
        assert "-b" in argv and "-y" in argv and "-I100" in argv
        assert "-s" not in argv
    second = engine.fit([9, 3, 4], True)
    assert first == second
    assert len(list(engine.cache.glob("*/attempt_*/command.json"))) == 2
    # Changing native outputs cannot silently alter a resumed analysis.
    (Path(first.run_directory) / "cafe.log").write_text("tampered")
    with pytest.raises(m.NativeCafeError, match="Cached"):
        engine.fit([9, 3, 4], True)


def test_restart_disagreement_is_not_convergence(tree_file, native_stub, tmp_path, monkeypatch):
    engine = m.NativeCafe(m.read_cafe_tree(tree_file, "<5>"), tmp_path / "work", str(native_stub), 2, 100, 1, 10)
    original = engine.command
    calls = []
    def disagree(directory, args):
        original(directory, args)
        calls.append(directory)
        (directory / "cafe.log").write_text(f"Completed 20 iterations\nFinal -lnL: {len(calls)}\n")
    monkeypatch.setattr(engine, "command", disagree)
    with pytest.raises(m.NativeCafeError, match="disagree"):
        engine.fit([9, 3, 4], False)


def test_failure_is_preserved_not_changed_to_p_one(tmp_path):
    class Engine:
        def fit(self, counts, alternative):
            if counts == [2]:
                raise m.NativeCafeError("replicate failed")
            return fit(12 if not alternative else 9, (.1, .2) if alternative else (.1,))
        def null_root(self, *args):
            return 4, "root"
        def simulate(self, *args):
            return [[2]], "sim"
    with pytest.raises(m.NativeCafeError, match="replicate failed"):
        m.compare_family(Engine(), [1], 1, tmp_path)
    assert (tmp_path / "observed.json").exists()
    assert not (tmp_path / "bootstrap.tsv").exists()


@pytest.mark.parametrize("record", [{}, [], {"directory": "../outside", "sha256": {}},
                                    {"directory": "attempt_test", "sha256": []}])
def test_corrupt_cache_manifest_fails_with_an_auditable_error(tree_file, native_stub, tmp_path, record):
    engine = m.NativeCafe(m.read_cafe_tree(tree_file, "<5>"), tmp_path / "work", str(native_stub), 2, 100, 1, 10)
    engine.fit([9, 3, 4], True)
    next(engine.cache.glob("*/complete.json")).write_text(json.dumps(record))
    with pytest.raises(m.NativeCafeError, match="Malformed"):
        engine.fit([9, 3, 4], True)


def test_changed_native_input_invalidates_saved_fit(tree_file, native_stub, tmp_path):
    engine = m.NativeCafe(m.read_cafe_tree(tree_file, "<5>"), tmp_path / "work", str(native_stub), 2, 100, 1, 10)
    fitted = engine.fit([9, 3, 4], True)
    (Path(fitted.run_directory) / "family.tsv").write_text("different counts")
    with pytest.raises(m.NativeCafeError, match="Cached"):
        engine.fit([9, 3, 4], True)


def test_competing_writer_does_not_replace_active_results(tmp_path):
    metadata = tmp_path / "metadata.json"
    metadata.write_text('{"status": "running"}')
    args = SimpleNamespace(output_dir=tmp_path)
    with m.output_lock(tmp_path):
        with pytest.raises(m.OutputBusyError):
            m.run(args)
    assert json.loads(metadata.read_text()) == {"status": "running"}
    with m.output_lock(tmp_path):
        pass  # Normal scope exit releases ownership.


def run_arguments(tree_file, tmp_path):
    args = SimpleNamespace(output_dir=tmp_path / 'out', asr_tree=tree_file,
        target_branch='A<1>', tree=None, counts=tmp_path / 'counts.tsv',
        changes=tmp_path / 'changes.tsv', family_ids=tmp_path / 'families.tsv',
        cafe='unused', fit_restarts=2, max_iterations=100, cores=1,
        timeout=10, error_model=None, bootstrap_replicates=1)
    labels = ['A<1>', 'B<2>', 'C<3>', '<4>', '<5>']
    m.write_table(args.counts, [dict(zip(['FamilyID', *labels], ['F', 9, 3, 4, 4, 4]))], ['FamilyID', *labels])
    m.write_table(args.changes, [dict(zip(['FamilyID', *labels], ['F', 5, -1, 0, 0, 0]))], ['FamilyID', *labels])
    args.family_ids.write_text('FamilyID\nF\n')
    return args


def test_invalid_family_header_records_failure(tree_file, tmp_path):
    args = run_arguments(tree_file, tmp_path)
    args.family_ids.write_text('wrong_header\nF\n')
    assert m.run(args) == 1
    metadata = json.loads((args.output_dir / 'metadata.json').read_text())
    assert metadata['status'] == 'failed' and 'FamilyID column' in metadata['error']


def test_failed_rerun_removes_prior_family_audit(tree_file, tmp_path, monkeypatch):
    args = run_arguments(tree_file, tmp_path)
    audit = args.output_dir / 'families' / m.hashlib.sha256(b'F').hexdigest()
    audit.mkdir(parents=True)
    for name in ('observed.json', 'generation.json', 'bootstrap.tsv'):
        (audit / name).write_text('previous successful analysis')
    monkeypatch.setattr(m, 'NativeCafe', lambda *a: SimpleNamespace(identity={}))
    def fail(*a):
        raise m.NativeCafeError('native fit failed')
    monkeypatch.setattr(m, 'compare_family', fail)
    assert m.run(args) == 1
    assert not list(audit.iterdir())
    _, rows = m.read_table(args.output_dir / 'family_lrt.tsv')
    assert rows[0]['status'] == 'failed' and rows[0]['p_value'] == 'NA'

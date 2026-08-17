import importlib
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SUPPORT_DIR = REPO_ROOT / "workflow" / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

import format_species_annotations as annotations  # noqa: E402
import format_species_download_runtime as download_runtime  # noqa: E402
import format_species_download_stage as download_stage  # noqa: E402
import format_species_inputs as format_species  # noqa: E402
import format_species_provider_inputs as provider_inputs  # noqa: E402
import format_species_provider_resolvers as provider_resolvers  # noqa: E402
import format_species_writers as writers  # noqa: E402


def test_format_species_inputs_delegates_writer_helpers():
    assert format_species.write_fasta_records_gzip is writers.write_fasta_records_gzip
    assert format_species.write_gff_gzip is writers.write_gff_gzip
    assert format_species.write_text_lines is writers.write_text_lines


def test_format_species_inputs_delegates_stage_and_provider_input_helpers():
    assert format_species.apply_download_input_dir is download_stage.apply_download_input_dir
    assert format_species.run_download_stage is download_stage.run_download_stage
    assert format_species.manifest_declared_providers is provider_inputs.manifest_declared_providers
    assert format_species.manifest_declared_species_keys is provider_inputs.manifest_declared_species_keys
    assert format_species.resolve_provider_inputs is provider_inputs.resolve_provider_inputs


def test_format_species_facade_preserves_the_legacy_public_import_surface():
    wildcard_modules = (
        "format_species_annotations",
        "format_species_cli",
        "format_species_common",
        "format_species_constants",
        "format_species_discovery",
        "format_species_download_runtime",
    )
    trailing_wildcard_modules = (
        "format_species_manifest",
        "format_species_provider_resolvers",
        "format_species_provider_urls",
        "format_species_summary",
        "format_species_taxonomy",
        "format_species_writers",
    )
    expected = {
        "sys": sys,
        "Path": Path,
        "SUPPORT_DIR": SUPPORT_DIR,
        "json": importlib.import_module("json"),
    }

    def add_wildcard_exports(module_name):
        module = importlib.import_module(module_name)
        names = getattr(module, "__all__", tuple(name for name in vars(module) if not name.startswith("_")))
        expected.update({name: getattr(module, name) for name in names})

    for module_name in wildcard_modules:
        add_wildcard_exports(module_name)
    expected.update(
        {
            "apply_download_input_dir": download_stage.apply_download_input_dir,
            "run_download_stage": download_stage.run_download_stage,
        }
    )
    add_wildcard_exports(trailing_wildcard_modules[0])
    expected.update(
        {
            "manifest_declared_providers": provider_inputs.manifest_declared_providers,
            "manifest_declared_species_keys": provider_inputs.manifest_declared_species_keys,
            "resolve_provider_inputs": provider_inputs.resolve_provider_inputs,
        }
    )
    for module_name in trailing_wildcard_modules[1:]:
        add_wildcard_exports(module_name)
    expected["main"] = format_species.main

    current_public_names = {name for name in vars(format_species) if not name.startswith("_")}
    assert current_public_names == set(expected)
    assert not hasattr(format_species, "__all__")
    for name, value in expected.items():
        if name == "SUPPORT_DIR":
            assert getattr(format_species, name) == value, name
        else:
            assert getattr(format_species, name) is value, name


def test_format_species_facade_stays_small_and_implementation_is_partitioned():
    facade_lines = (SUPPORT_DIR / "format_species_inputs.py").read_text(encoding="utf-8").splitlines()
    expected_modules = {
        "format_species_annotations.py",
        "format_species_cli.py",
        "format_species_common.py",
        "format_species_constants.py",
        "format_species_discovery.py",
        "format_species_download_runtime.py",
        "format_species_provider_resolvers.py",
        "format_species_run_summary.py",
        "format_species_summary.py",
        "format_species_taxonomy.py",
    }

    assert len(facade_lines) < 400
    assert expected_modules <= {path.name for path in SUPPORT_DIR.glob("format_species_*.py")}


def test_large_format_species_facades_delegate_to_bounded_components():
    facade_modules = {
        annotations: "format_species_annotation",
        download_runtime: "format_species_download",
        provider_resolvers: "format_species_providers",
    }
    for facade, package_name in facade_modules.items():
        facade_path = Path(facade.__file__)
        component_dir = SUPPORT_DIR / package_name
        component_paths = sorted(path for path in component_dir.glob("*.py") if path.name != "__init__.py")

        assert len(facade_path.read_text(encoding="utf-8").splitlines()) < 300
        assert component_paths
        assert all(len(path.read_text(encoding="utf-8").splitlines()) < 900 for path in component_paths)

    assert annotations.derive_cds_records_from_gbff.__module__ == "format_species_annotation.genbank"
    assert download_runtime.download_url_to_file.__module__ == "format_species_download.locking"
    assert provider_resolvers.resolve_ncbi_download_urls_from_id.__module__ == "format_species_providers.ncbi"
    assert provider_resolvers.resolve_fernbase_download_urls_from_id.__module__ == (
        "format_species_providers.fernbase_oryza"
    )

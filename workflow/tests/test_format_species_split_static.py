import importlib
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SUPPORT_DIR = REPO_ROOT / "workflow" / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

import format_species_download_stage as download_stage  # noqa: E402
import format_species_inputs as format_species  # noqa: E402
import format_species_provider_inputs as provider_inputs  # noqa: E402


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

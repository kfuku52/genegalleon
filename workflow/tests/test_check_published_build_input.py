from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "container"
    / "scripts"
    / "check_published_build_input.py"
)


def load_module():
    spec = spec_from_file_location("check_published_build_input", SCRIPT)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def image(runtime_input="a" * 64, revision="b" * 40, build_target="runtime"):
    return {
        "os": "linux",
        "architecture": "amd64",
        "config": {
            "Labels": {
                "io.genegalleon.runtime-input": runtime_input,
                "io.genegalleon.build-target": build_target,
                "io.genegalleon.security-refresh-epoch": "2026-08-31",
                "org.opencontainers.image.revision": revision,
            }
        },
    }


def test_single_platform_image_matches_only_its_declared_platform():
    module = load_module()
    payload = image()

    assert module.compare_build_inputs(
        payload,
        {"linux/amd64": "a" * 64},
        label="io.genegalleon.runtime-input",
    ) == []
    assert module.compare_build_inputs(
        payload,
        {"linux/arm64": "a" * 64},
        label="io.genegalleon.runtime-input",
    ) == ["Published image is missing platform linux/arm64."]


def test_single_platform_image_exposes_its_common_revision():
    module = load_module()

    assert module.common_revision(image(), ["linux/amd64"]) == ("b" * 40, [])


def test_multi_platform_mapping_remains_unchanged():
    module = load_module()
    payload = {"linux/amd64": image()}

    assert module.platform_mapping(payload) is payload


def test_single_platform_exact_contract_labels_are_compared():
    module = load_module()

    assert module.compare_label_values(
        image(),
        {"linux/amd64": "runtime"},
        label="io.genegalleon.build-target",
    ) == []
    assert module.compare_label_values(
        image(build_target="development"),
        {"linux/amd64": "runtime"},
        label="io.genegalleon.build-target",
    ) == [
        "Published image linux/amd64 io.genegalleon.build-target changed: "
        "published=development expected=runtime."
    ]

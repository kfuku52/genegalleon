from shell_static_helpers import REPO_ROOT, read_text


def test_amd64_container_pins_trinity_compatible_salmon():
    requirements = read_text(REPO_ROOT / "container" / "env" / "base.required.txt")
    commands = read_text(REPO_ROOT / "container" / "spec" / "required_commands.tsv")

    assert "salmon=1.11.4" in requirements.splitlines()
    assert "base\tsalmon" in commands.splitlines()


def test_container_builds_generate_salmon_locale():
    dockerfile = read_text(REPO_ROOT / "container" / "Dockerfile")
    definition = read_text(REPO_ROOT / "container" / "apptainer_local_build.def.template")

    for build_file in (dockerfile, definition):
        assert "locales" in build_file
        assert "locale-gen en_US.UTF-8" in build_file


def test_runtime_validation_exercises_trinity_salmon_interface():
    validation = read_text(REPO_ROOT / "container" / "scripts" / "validate_runtime.sh")

    assert "check_trinity_salmon_compatibility" in validation
    assert "Trinity requires Salmon 1.x" in validation
    assert "--minAssignedFrags 1" in validation
    assert "--validateMappings" in validation

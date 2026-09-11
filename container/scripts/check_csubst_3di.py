#!/usr/bin/env python3
"""Check the installed CSUBST 3Di contract without downloading model weights."""

import importlib
import json
from importlib.metadata import metadata, version

from packaging.requirements import Requirement


def main():
    from csubst import scan_endpoint
    from csubst.recoding_config import DEFAULT_SA_BACKEND
    from peft import LoraConfig, TaskType, get_peft_model
    from transformers import EsmConfig, EsmForTokenClassification, EsmTokenizer, T5EncoderModel, T5Tokenizer

    scan_endpoint.validate_options({
        "subcommand": "scan", "nonsyn_recode": "3di20", "sa_asr_mode": "direct",
        "sa_iqtree_model": "GTR", "scan_observation": "joint", "scan_rate_exposure": "endpoint",
        "scan_rate_length": "raw", "scan_rate_event_mode": "posterior_sum", "scan_pvalue_calibration": "none",
    })
    # Keep the runtime check tied to the owning package's optional dependencies.
    dependencies = {}
    for value in metadata("csubst").get_all("Requires-Dist", []):
        requirement = Requirement(value)
        if requirement.marker and requirement.marker.evaluate({"extra": "3di"}):
            installed = version(requirement.name)
            if installed not in requirement.specifier:
                raise RuntimeError(f"CSUBST 3Di requires {requirement}; found {installed}")
            dependencies[requirement.name] = installed
    for name in ("torch", "huggingface_hub", "sentencepiece", "google.protobuf"):
        importlib.import_module(name)
    assert all((LoraConfig, TaskType, get_peft_model, EsmConfig, EsmForTokenClassification,
                EsmTokenizer, T5EncoderModel, T5Tokenizer))
    print(json.dumps({"csubst": version("csubst"), "default_backend": DEFAULT_SA_BACKEND,
                      "dependencies": dependencies, "model_download": False}, sort_keys=True))


if __name__ == "__main__":
    main()

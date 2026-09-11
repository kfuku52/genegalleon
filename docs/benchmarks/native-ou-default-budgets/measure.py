"""Sequential fresh-process wall/CPU/RSS timing with a one-hour safety cutoff."""

import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "existing_measure", ROOT.parent / "native-ou-1000tips-100shifts-aicc" / "benchmark.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
label, traits, pool, refits, screens = sys.argv[1:]
output = ROOT / "runs" / label
output.mkdir(parents=True, exist_ok=False)
command = [
    sys.executable,
    str(ROOT / "run_search.py"),
    str(ROOT / "data" / f"1000tips-{traits}traits"),
    str(output / "model.json"),
    "--candidate-pool",
    pool,
    "--refit-budget",
    refits,
    "--screening-budget",
    screens,
]
result = module.measure(command, str(output / "run.log"), 3600)
(output / "measurement.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"label": label, **result}), flush=True)

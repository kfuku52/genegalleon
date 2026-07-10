#!/usr/bin/env python3

import runpy
from pathlib import Path

if __name__ == "__main__":
    target = Path(__file__).with_name("summarize_gg_input_generation_runs.py")
    runpy.run_path(str(target), run_name="__main__")

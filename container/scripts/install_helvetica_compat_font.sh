#!/usr/bin/env bash
set -euo pipefail

env_name=${1:-base}

log() {
  echo "[install_helvetica_compat_font] $*"
}

main() {
  micromamba run -n "${env_name}" python - <<'PY'
from pathlib import Path

import matplotlib
from fontTools.ttLib import TTFont

ttf_dir = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
src = ttf_dir / "DejaVuSans.ttf"
dst = ttf_dir / "Helvetica.ttf"

if not src.exists():
    raise SystemExit(f"Source font not found: {src}")

font = TTFont(str(src))
for rec in font["name"].names:
    if rec.nameID in (1, 4, 16):
        text = "Helvetica"
    elif rec.nameID in (2, 17):
        text = "Regular"
    elif rec.nameID == 6:
        text = "Helvetica-Regular"
    else:
        continue

    if rec.platformID in (0, 3):
        rec.string = text.encode("utf_16_be")
    else:
        rec.string = text.encode("latin-1", errors="ignore")

font.save(str(dst))
print(f"Installed {dst}")
PY

  rm -rf -- "/root/.cache/matplotlib" "/home/.cache/matplotlib" || true

  micromamba run -n "${env_name}" python - <<'PY'
import matplotlib
from matplotlib import font_manager

font_manager.fontManager = font_manager._load_fontmanager(try_read_cache=False)
resolved = font_manager.findfont("Helvetica")
print(f"Matplotlib Helvetica path: {resolved}")

if not resolved.endswith("Helvetica.ttf"):
    raise SystemExit(f"Failed to install Helvetica compatibility font. Resolved path: {resolved}")
PY

  log "Helvetica compatibility font installed for Matplotlib."
}

main "$@"

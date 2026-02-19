#!/usr/bin/env bash
set -euo pipefail

arch=${1:-unknown}
note_file=/opt/pg/logs/BUILD_NOTES.txt

cat > "${note_file}" <<EOF
genegalleon multi-arch build notes
==========================
Target arch: ${arch}

This image is a reproducible scaffold for genegalleon and is intentionally conservative.
The original genegalleon.sif was assembled interactively and includes tools/databases that may be:
1) not available on all architectures, or
2) too large to redistribute in a public Dockerfile build.

Expected manual post-build steps:
- Populate /usr/local/db/Pfam_LE
- Populate /usr/local/db/uniprot_sprot.pep (and DIAMOND DB if required)
- Populate /usr/local/db/jaspar
- Verify Notung jar exists at /usr/local/bin/Notung.jar
- Validate workflow/gg_test_cmd.sh in each architecture image

Optional package installation logs:
- /opt/pg/logs/failed_optional_base.txt

Runtime validation report:
- /opt/pg/logs/runtime_validation_${arch}.tsv
EOF

echo "[post_install_notes] Wrote ${note_file}"

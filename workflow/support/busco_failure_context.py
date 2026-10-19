"""Record a failed BUSCO attempt without inferring biological completeness."""
import argparse
import json
from pathlib import Path


def capture(root, mode, lineage, exit_code):
    root = Path(root)
    paths = [root / 'input.fasta']
    for suffix in ('fas', 'codon.fas', 'modified.fas'):
        paths.extend((root / 'busco_tmp' / ('run_' + lineage) / 'metaeuk_output' / stage /
                      ('input.fasta.' + suffix)) for stage in ('initial_results', 'rerun_results'))
    files = []
    for path in paths:
        if path.exists():
            files.append(dict(path=str(path.relative_to(root)), size=path.stat().st_size,
                              symlink=path.is_symlink()))
    value = dict(schema_version=1, status='failed', mode=mode, lineage=lineage,
                 exit_code=exit_code, files=files)
    # Absence of a completed full/short table is never a zero-percent result.
    (root / 'failure-context.json').write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--mode', choices=('transcriptome', 'proteins'), required=True)
    parser.add_argument('--lineage', required=True)
    parser.add_argument('--exit-code', type=int, required=True)
    args = parser.parse_args()
    capture(args.root, args.mode, args.lineage, args.exit_code)


if __name__ == '__main__':
    main()

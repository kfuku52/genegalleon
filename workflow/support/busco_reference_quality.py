"""Pure BUSCO quality, tree-distance and identifier contracts.

These helpers do not require the optional genome-rescue/synteny toolchain.
"""
import math
import re
from collections import defaultdict
from pathlib import Path

COMPARABLE_QUALITY = ("lineage", "version", "mode", "lineage_date", "markers")


def busco_quality(path):
    text = Path(path).read_text()
    complete = re.search(r"C:([\d.]+)%", text)
    lineage = re.search(r"lineage dataset is:\s*(\S+)", text)
    version = re.search(r"BUSCO version is:\s*(\S+)", text)
    mode = re.search(r"BUSCO was run in mode:\s*(\S+)", text)
    markers = re.search(r"\bn:\s*(\d+)", text)
    date = re.search(r"Creation date:\s*([^,\s)]+)", text)
    if not all((complete, lineage, version, mode, markers)):
        raise ValueError(f"BUSCO summary lacks completeness/lineage/version/mode/marker count: {path}")
    value = float(complete[1])
    if not math.isfinite(value) or not 0 <= value <= 100:
        raise ValueError(f"Invalid BUSCO completeness: {path}")
    if int(markers[1]) < 1:
        raise ValueError(f"Invalid BUSCO marker count: {path}")
    return {"complete_pct": value, "lineage": lineage[1], "version": version[1], "mode": mode[1],
            "lineage_date": date[1] if date else None, "markers": int(markers[1])}



def patristic_distances(tree):
    """All tip distances in O(N^2), without repeated whole-tree LCA scans."""
    graph = defaultdict(list)
    for parent in tree.find_clades():
        for child in parent.clades:
            length = child.branch_length or 0.0
            graph[parent].append((child, length))
            graph[child].append((parent, length))
    leaves = {tip: tip.name for tip in tree.get_terminals()}
    result = {}
    for tip, name in leaves.items():
        row, stack = {}, [(tip, None, 0.0)]
        while stack:
            node, previous, distance = stack.pop()
            if node in leaves:
                row[leaves[node]] = distance
            stack.extend((child, node, distance + length) for child, length in graph[node] if child is not previous)
        result[name] = row
    return result



def safe_token(value, label):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError(f"{label} must be a safe identifier: {value!r}")
    return value


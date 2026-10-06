"""Record and replay the gg_gene_evolution tree renderer's argument contract.

Older results have parameter/input provenance rather than an argument snapshot;
their panel list is reconstructed explicitly, with the source recorded in the
focused output. Missing optional measurements remain unavailable to treevis.
"""

import argparse
import gzip
import json
import re
from pathlib import Path


def record(output, family_root, arguments, species_parser):
    roots = {str(Path(family_root).resolve()), str(Path(family_root).absolute())}
    for root in sorted(roots, key=len, reverse=True):
        arguments = [arg.replace(root + '/', '{family_root}/') for arg in arguments]
    payload = dict(
        schema_version=1,
        arguments=arguments,
        species_label_parser=species_parser,
    )
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(json.dumps(payload, indent=2) + "\n")


def fasta_ids(raw, name):
    if name.endswith(".gz"):
        raw = gzip.decompress(raw)
    return {line[1:].split()[0] for line in raw.decode().splitlines() if line.startswith(">")}


def replay(store, family, rows, destination, sources):
    """Return original settings with all available family inputs materialized."""
    import hashlib

    def read(subdir, name):
        try:
            with store.open_binary(subdir, name) as h:
                raw = h.read()
        except FileNotFoundError:
            return None
        sources[subdir + "/" + name] = hashlib.sha256(raw).hexdigest()
        return raw

    snapshot = read("artifact_provenance", family + ".tree_plot.args.json")
    provenance = read("artifact_provenance", family + ".tree_plot.json")
    metadata = json.loads(provenance) if provenance else {}
    parameters = metadata.get("parameters", {})
    paths = {}
    # Recorded input paths select the actual analysis/pruned alignment, not a
    # merely equal-sized family alignment. Legacy paths are read by the store.
    for item in metadata.get("inputs", []):
        if item.get("scope") == "logical":
            paths[item["label"]] = item["path"]
    optional = {
        "synteny": f"synteny/{family}_synteny.tsv",
        "domain": f"rpsblast/{family}_rpsblast.tsv",
        "trimmed": paths.get("input_3", f"clipkit/{family}_cds.clipkit.fa.gz"),
        "untrimmed": paths.get("input_4", f"mafft/{family}_cds.aln.fa.gz"),
        "dated": f"dated_tree/{family}_dated.nwk",
        "fimo": f"fimo/{family}_fimo.tsv",
        "meme": f"meme/{family}_meme.xml",
        "promoter": f"promoter_fasta/{family}_promoter.fa.gz",
    }
    raw_files = {}

    def materialize(logical):
        if logical in raw_files:
            return str(destination / logical)
        if logical.startswith("/") or ".." in Path(logical).parts:
            raise ValueError("Unsafe family renderer input: " + logical)
        subdir, name = logical.split("/", 1)
        raw = read(subdir, name)
        path = destination / logical
        if raw is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(raw)
            raw_files[logical] = raw
        return str(path)

    for logical in optional.values():
        materialize(logical)
    for logical in paths.values():
        materialize(logical)
    tips = {r["node_name"] for r in rows if r.get("child1") == r.get("child2") == "-999"}
    for kind, candidates in (
        (
            "trimmed",
            [
                optional["trimmed"],
                f"orthogroup_extraction_fasta/{family}_orthogroup_extraction.fa.gz",
                f"maxalign/{family}_cds.maxalign.fa.gz",
                f"mafft/{family}_cds.aln.fa.gz",
            ],
        ),
        (
            "untrimmed",
            [
                optional["untrimmed"],
                f"orthogroup_extraction_fasta/{family}_orthogroup_extraction.fa.gz",
                f"mafft/{family}_cds.aln.fa.gz",
            ],
        ),
    ):
        for logical in dict.fromkeys(candidates):
            materialize(logical)
            if logical in raw_files and tips <= fasta_ids(raw_files[logical], logical):
                optional[kind] = logical
                break
        else:
            optional[kind] = ""
    if snapshot:
        spec = json.loads(snapshot)
        if spec.get("schema_version") != 1:
            raise ValueError("Unsupported gene-tree renderer argument contract")
        arguments = []
        for arg in spec["arguments"]:
            for logical in re.findall(r"\{family_root\}/([^,|\s]+)", arg):
                materialize(logical)
            arguments.append(arg.replace("{family_root}", str(destination)))
        parser = spec["species_label_parser"]
        source = "recorded_gg_gene_evolution_arguments"
    else:
        if not provenance:
            raise ValueError("gg_gene_evolution renderer settings unavailable for " + family)

        def p(key, default):
            return str(parameters.get(key, default))

        def file(key):
            return str(destination / optional[key]) if optional[key] else ""

        panels = [
            f"tree,{p('branch_length', 'bl_rooted')},{p('support_value_resolved', 'no')},{p('branch_color', 'no')},L",
            f"heatmap,{p('heatmap_transform', 'no')},abs,_,expression_",
            "pointplot,no,rel,_,expression_",
            f"cluster_membership,{p('max_intergenic_dist', 100000)}",
            f"synteny,{file('synteny')},{p('synteny_window', 5)}",
            "tiplabel",
        ]
        if any(row.get("query_marker") for row in rows) and p("query_marker", 1) == "1":
            panels += ["categorical,query_marker,Query,-"]
        cds = file("untrimmed") if p("sequence_similarity_mode", "cds") == "cds" else ""
        prefixes = {
            str(row.get("ortholog", "")).rsplit("_", 1)[0] + "_"
            for row in rows
            if row.get("ortholog") and row.get("ortholog") not in {"NA", "nan"}
        }
        prefix = p("clade_ortholog_prefix", "")
        if not prefix:
            prefix = sorted(prefixes)[0] if len(prefixes) == 1 else "HGT_UNRESOLVED_"
        if p("clade_ortholog", 1) != "1":
            prefix = ""
        panels += [
            "localization",
            "transmembrane_domain",
            "intron_number",
            f"gene_structure,compressed,23,{cds}",
            f"domain,{file('domain')}",
            f"alignment,{file('trimmed')},{file('untrimmed')}",
            f"fimo,{p('promoter_bp', 2000)},{p('fimo_qvalue', 0.05)}",
            f"meme,{file('meme')}",
            f"ortholog,{prefix},{file('dated')}",
        ]
        if p("synteny_similarity", 1) == "1":
            panels += [
                f"synteny_similarity,{file('synteny')},{p('synteny_search_window', 20)},{p('synteny_similarity_width_mm', 15)}"
            ]
        if p("cis_similarity", 1) == "1":
            panels += [
                f"cis_similarity,{file('fimo')},{file('promoter')},{p('fimo_qvalue', 0.05)},{p('cis_similarity_width_mm', 15)}"
            ]
        if p("sequence_similarity", 1) == "1":
            panels += [
                f"sequence_similarity,{file('trimmed')},{p('sequence_similarity_mode', 'cds')},{p('sequence_similarity_width_mm', 15)},{p('sequence_similarity_genetic_code', 1)}"
            ]
        arguments = [f"--panel{i}={panel}" for i, panel in enumerate(panels, 1)]
        arguments += [
            "--panel_widths_mm=tree:60",
            "--show_branch_id=yes",
            "--species_color_table=PLACEHOLDER",
            f"--max_delta_intron_present={p('retrotransposition_delta_intron', -0.5)}",
            f"--event_method={p('event_method', 'auto')}",
            f"--pie_chart_value_transformation={p('pie_chart_value_transformation', 'identity')}",
        ]
        for key, default in [
            ("long_branch_display", "auto"),
            ("long_branch_ref_quantile", 0.95),
            ("long_branch_detect_ratio", 5),
            ("long_branch_cap_ratio", 2.5),
            ("long_branch_tail_shrink", 0.02),
            ("long_branch_max_fraction", 0.1),
        ]:
            arguments.append(f"--{key}={p(key, default)}")
        # Protein convergence is optional and no CB result is invented.
        arguments += [
            f"--protein_convergence=100,100,yes,3-{p('csubst_max_arity', 10)},"
            f"{destination}/csubst_cb/{family}_cb_ARITY.tsv,{p('csubst_cutoff_stat', '')}"
        ]
        parser = p("species_label_parser", "taxonomic")
        source = "saved_gg_gene_evolution_parameter_provenance"
    arguments = [arg for arg in arguments if not arg.startswith("--stat_branch=")]
    panels = [arg for arg in arguments if re.match(r"--panel\d+=", arg)]
    last = max([int(arg.split("=")[0][7:]) for arg in panels], default=0)
    arguments.append(f"--panel{last + 1}=categorical,hgtfocus_tip_status,Scaffold-supported descendants,-")
    return dict(
        arguments=arguments,
        species_label_parser=parser,
        settings_source=source,
        optional_input_availability={key: value in raw_files for key, value in optional.items()},
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--family-root", required=True)
    parser.add_argument("--species-parser", required=True)
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    record(args.output, args.family_root, [x for x in args.arguments if x != "--"], args.species_parser)


if __name__ == "__main__":
    main()

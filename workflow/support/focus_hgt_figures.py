"""Three source-backed, single-page summaries of supported category-1 HGT."""

import csv
import gzip
import io
import math
from collections import Counter

import numpy as np
from gene_family_output_store import GeneFamilyOutputStore, read_only_observation

ORANGE = "#b34d00"
BLUE = "#2b6ca3"


def read(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write(path, rows, fields=None):
    with path.open("w", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=fields or list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def filtering_counts(source_events, selected, audit_path="", prefilter_selected=None, pfam_selected=None,
                     direction_selected=None, support_filter_enabled=False):
    def unique(rows):
        result = {row['event_id']: row for row in rows}
        if len(result) != len(rows) or '' in result:
            raise ValueError('Duplicate or empty modeled event ID in filtering cohort')
        return result

    def subset(rows, parent):
        for row in rows:
            original = parent.get(row['event_id'])
            if original is None:
                raise ValueError('Focused input cohort is not a subset of its source events')
            for field in ('orthogroup', 'gene_tree_branch_id', 'branch_id', 'gene_tree_node', 'node_name',
                          'generax_donor_node', 'generax_recipient_node', 'generax_transfer', 'event_index'):
                if field in row and field in original and row[field] != original[field]:
                    raise ValueError('Filtering event identity disagrees: ' + field)
            for aliases in [('gene_tree_branch_id', 'branch_id'), ('gene_tree_node', 'node_name')]:
                identities = [{str(record[k]) for k in aliases if record.get(k) not in (None, '')}
                              for record in (row, original)]
                if any(len(values) > 1 for values in identities) or all(identities) and identities[0] != identities[1]:
                    raise ValueError('Filtering event identity disagrees: ' + '/'.join(aliases))

    source_by_id = unique(source_events)
    unique(selected)
    subset(selected, source_by_id)
    if prefilter_selected is not None and pfam_selected is not None:
        raise ValueError('Supply only one Pfam/trait filtering order')
    if pfam_selected is not None:
        pfam_by_id = unique(pfam_selected)
        subset(pfam_selected, source_by_id)
        subset(selected, pfam_by_id)
    if prefilter_selected is not None:
        prefilter_by_id = unique(prefilter_selected)
        subset(prefilter_selected, source_by_id)
        subset(selected, prefilter_by_id)
    if direction_selected is not None:
        if prefilter_selected is not None:
            raise ValueError('Species-branch direction filtering requires Pfam-before-trait order')
        direction_by_id = unique(direction_selected)
        subset(direction_selected, pfam_by_id if pfam_selected is not None else source_by_id)
        subset(selected, direction_by_id)
    stages = []
    if audit_path:
        audited = read(audit_path)
        ids = [row["event_id"] for row in audited]
        if len(ids) != len(set(ids)):
            raise ValueError("Duplicate modeled event in filtering audit")
        if not support_filter_enabled:
            audited_by_id = unique(audited)
            subset(source_events, audited_by_id)
            for row in source_events:
                original = audited_by_id[row['event_id']]
                for field, expected in (('mapping_status', 'matched'),
                                        ('species_tree_mapping_status', 'matched_external_tree')):
                    if field in original and original[field] != expected:
                        raise ValueError('Focused cohort contains an unresolved event/species mapping')
            stages.append(("All modeled transfers", audited))
        else:
            directional = [
                row
                for row in audited
                if row.get("donor_classification") == "outside" and row.get("recipient_classification") == "insect"
            ]
            accepted = [row for row in audited if row.get("status") == "accepted"]
            for row in accepted:
                value = float(row.get("support_used") or row.get("support_generax_ufboot") or "nan")
                support_source = row.get("support_source", "")
                if (
                    not math.isfinite(value)
                    or not 90 <= value <= 100
                    or not ("support_generax_ufboot" in support_source or "raw_unrooted_split_verified" in support_source)
                ):
                    raise ValueError("Accepted filtering-audit event lacks verified inclusive UFBoot >=90")
            accepted_ids = {row["event_id"] for row in accepted}
            if not {row["event_id"] for row in source_events} <= accepted_ids:
                raise ValueError("Focused input cohort is not a subset of accepted filtering-audit events")
            subset(source_events, {r['event_id']: r for r in accepted})
            if not accepted_ids <= {row["event_id"] for row in directional}:
                raise ValueError("Accepted filtering-audit event has unsupported transfer direction")
            stages.append(("All modeled transfers", audited))
            # Historical direction decisions stay in the source audit. Focused
            # figures show taxonomy only at the final, combined selection stage.
            if direction_selected is None:
                stages.append(("Non-Insecta to Insecta", directional))
            stages.append(("Matched gene-tree UFB >=90", accepted))
    stages.append(("Input supported-event cohort", source_events))
    if pfam_selected is not None:
        stages.append(("Event-gene pair Pfam filter", pfam_selected))
        stages.append(("Non-Arthropoda donor & category = 1 recipient"
                       if direction_selected is not None else "Category = 1 recipients", selected))
    else:
        # Compatibility for callers reproducing the earlier trait-first figure.
        stages.append(("Non-Arthropoda donor & category = 1 recipient"
                       if direction_selected is not None else "Category = 1 recipients",
                       selected if prefilter_selected is None else prefilter_selected))
        if prefilter_selected is not None:
            stages.append(("Event-gene pair Pfam filter", selected))
    return [
        dict(stage=label, event_count=len(rows), orthogroup_count=len({r["orthogroup"] for r in rows}))
        for label, rows in stages
    ]


def available_label(value):
    return (
        value.strip() if value and value.strip().lower() not in {"na", "nan", "none", "annotation unavailable"} else ""
    )


def best_hit_product(row):
    name = available_label(row.get('swissprot_best_hit_protein_name', ''))
    if name:
        return name
    basis = row.get('best_available_product_label_basis', '').lower().replace('-', '').replace('_', '').replace(' ', '')
    if 'besthit' in basis and 'gff' not in basis:
        return available_label(row.get('best_available_product_label', ''))
    return ''


def product_labels(families, links):

    labels = {}
    for family in families:
        candidates = sorted(
            [r for r in links if r["orthogroup"] == family and r["side"] == "recipient"
             and str(r.get('eligible_for_context', 'True')).lower() in {'true', '1'}], key=lambda r: r["gene_id"]
        )
        known = [r for r in candidates if best_hit_product(r)]
        row = known[0] if known else {}
        labels[family] = dict(
            protein_product=best_hit_product(row) or "Annotation unavailable",
            annotation_gene_id=row.get("gene_id", ""),
            annotation_basis='SwissProt_best_hit_prediction'
            if available_label(row.get('swissprot_best_hit_protein_name', ''))
            else available_label(row.get("best_available_product_label_basis", "")),
            all_recipient_product_labels="; ".join(sorted({best_hit_product(r) for r in candidates if best_hit_product(r)})),
        )
    return labels


def export_filtering_flow(directory, source_events, selected, trait, filter_audit='',
                          prefilter_selected=None, pfam_selected=None, direction_selected=None,
                          support_filter_enabled=False, analyzed_orthogroups=None):
    """Render shared Pfam followed by one combined taxonomy/trait stage."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    directory.mkdir(parents=True, exist_ok=True)
    def save(fig, name, title, subtitle, note):
        fig.text(0.035, 0.965, title, fontsize=17, weight='bold', va='top')
        fig.text(0.035, 0.916, subtitle, fontsize=9, color='#556975', va='top')
        fig.text(0.035, 0.045, note, fontsize=8, color='#556975', va='bottom')
        fig.savefig(directory / (name + '.pdf'))
        fig.savefig(directory / (name + '.png'), dpi=140)
        plt.close(fig)
    counts = filtering_counts(source_events, selected, filter_audit,
                              prefilter_selected=prefilter_selected, pfam_selected=pfam_selected,
                              direction_selected=direction_selected, support_filter_enabled=support_filter_enabled)
    coverage_values = {float(r['pfam_min_shared_query_coverage']) for r in (pfam_selected or [])
                       if r.get('pfam_min_shared_query_coverage') not in (None, '')}
    if len(coverage_values) > 1:
        raise ValueError('Mixed shared-Pfam coverage thresholds in filtering cohort')
    if coverage_values:
        minimum = coverage_values.pop()
        for row in counts:
            if row['stage'] == 'Event-gene pair Pfam filter':
                row['stage'] = f'Event-pair Pfam (>={100 * minimum:g}% each query)'
    trait_index = -1 if prefilter_selected is None else -2
    if filter_audit:
        for row in counts:
            if row['stage'] == 'Input supported-event cohort':
                row['stage'] = 'Bilateral scaffold background'
    counts[trait_index]["stage"] = ("Non-Arthropoda donor & " + trait + " = 1 recipient"
                                   if direction_selected is not None else trait + " = 1 recipients")
    write(directory / "filtering_flow_audit.tsv", counts)
    displayed = [dict(row) for row in counts if row['stage'] != 'Matched event and species branches']
    first_step = 1
    if analyzed_orthogroups is not None:
        families = list(analyzed_orthogroups)
        if len(families) != len(set(families)) or any(not family for family in families):
            raise ValueError('Duplicate or empty analyzed orthogroup ID')
        if not {row['orthogroup'] for row in source_events + selected} <= set(families):
            raise ValueError('Filtering cohort is outside analyzed orthogroups')
        if any(row['orthogroup_count'] > len(families) for row in counts):
            raise ValueError('Filtering orthogroup count exceeds analyzed orthogroups')
        displayed.insert(0, dict(stage='All analyzed orthogroups', event_count='NA', orthogroup_count=len(families)))
        first_step = 0
    displayed = [dict(step=f'{i + first_step:02d}', **row) for i, row in enumerate(displayed)]
    write(directory / "filtering_flow.tsv", displayed)
    fig, ax = plt.subplots(figsize=(13, 8))
    ax.set_axis_off()
    fig.subplots_adjust(top=0.80, bottom=0.14)
    for i, row in enumerate(displayed):
        y = 0.91 - i * 0.80 / max(1, len(displayed) - 1)
        color = ORANGE if i == len(displayed) - 1 else BLUE
        height = min(0.14, 0.72 / max(1, len(displayed) - 1))
        ax.add_patch(Rectangle((0.03, y - height / 2), 0.94, height, facecolor="#f0f4f7"))
        ax.text(0.06, y, row['step'], fontsize=16, color=color, va="center", weight="bold")
        ax.text(0.16, y, row["stage"], fontsize=13, va="center")
        event_label = f"{row['event_count']:,}" if isinstance(row['event_count'], int) else row['event_count']
        ax.text(0.76, y, event_label, ha="right", va="center", fontsize=20, color=color, weight="bold")
        ax.text(0.92, y, f"{row['orthogroup_count']:,}", ha="right", va="center", fontsize=15)
    ax.text(0.76, 1.09, "Events", ha="right", color="#556975")
    ax.text(0.92, 1.09, "Orthogroups", ha="right", color="#556975")
    save(
        fig,
        "filtering_flow",
        "From modeled transfers to focused candidates",
        "Distinct event IDs and orthogroups; duplication and repeated per-tip reports are not new modeled events."
        + ("\nStep 00 counts existing gene-tree summaries, including zero-transfer OGs; its event count is NA."
           if analyzed_orthogroups is not None else ""),
        ("No UFB threshold applied; measured support and missingness remain annotations.\n"
         if not support_filter_enabled else
         "UFB/scaffold/Pfam counts retain the previously verified input cohort; broader upstream support counts are not inferred.\n"
         if direction_selected is not None and filter_audit else
         "Upstream direction/support counts are shown only when an explicit event-level filtering audit is supplied.\n") +
        "UFB = Ultrafast bootstrap (gene-tree split support); category-1 internal recipient branches require all observed descendant tips = 1."
        + ("\nEvent/species correspondence checks remain in the audit."
           if not support_filter_enabled and filter_audit else ""),
    )
    return counts


def export_figures(directory, source_events, selected, links, tree, values, family_root, trait, filter_audit="",
                   context_annotations='', prefilter_selected=None, pfam_selected=None, direction_selected=None,
                   support_filter_enabled=False, analyzed_orthogroups=None):
    import hashlib
    import textwrap

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Rectangle

    directory.mkdir(parents=True)

    def save(fig, name, title, subtitle, note):
        fig.text(0.035, 0.965, title, fontsize=17, weight="bold", va="top")
        fig.text(0.035, 0.916, subtitle, fontsize=9, color="#556975", va="top")
        fig.text(0.035, 0.045, note, fontsize=8, color="#556975", va="bottom")
        fig.savefig(directory / (name + ".pdf"))
        fig.savefig(directory / (name + ".png"), dpi=140)
        plt.close(fig)

    counts = export_filtering_flow(directory, source_events, selected, trait, filter_audit,
                                   prefilter_selected=prefilter_selected, pfam_selected=pfam_selected,
                                   direction_selected=direction_selected, support_filter_enabled=support_filter_enabled,
                                   analyzed_orthogroups=analyzed_orthogroups)
    families = sorted({r["orthogroup"] for r in selected})
    species = [tip.name for tip in tree.get_terminals()]
    membership = Counter()
    sources = {}
    missing = set()
    from focus_hgt_gene_trees import validated_event_links
    annotation_links = [dict(r) for r in validated_event_links(selected, links)]
    from focus_hgt_context_annotations import ContextAnnotations, available
    annotations = ContextAnnotations(context_annotations)
    with read_only_observation():
        store = GeneFamilyOutputStore(family_root)
        with store.read_snapshot():
            for family in families:
                name = family + "_stat.branch.tsv"
                try:
                    with store.open_binary("stat_branch", name) as h:
                        raw = h.read()
                except FileNotFoundError:
                    missing.add(family)
                    continue
                sources["stat_branch/" + name] = hashlib.sha256(raw).hexdigest()
                rows = list(csv.DictReader(io.StringIO(raw.decode()), delimiter="\t"))
                tip_annotations = {row["node_name"]: row for row in rows if row.get('child1') == row.get('child2') == '-999'}
                if len(tip_annotations) != sum(row.get('child1') == row.get('child2') == '-999' for row in rows):
                    raise ValueError('Duplicate gene-tree tip in distribution annotations')
                for link in annotation_links:
                    if link['orthogroup'] != family:
                        continue
                    leaf = tip_annotations.get(link['gene_id'])
                    if leaf is None:
                        raise ValueError('Distribution annotation gene is absent from its family tree')
                    annotation = annotations.get(link['gene_id'], family, leaf)
                    if available(link.get('besthit_accession')) and 'sprot_best' in leaf \
                            and available(link['besthit_accession']) != available(leaf['sprot_best']):
                        raise ValueError('Distribution annotation best hit disagrees with the exact family leaf')
                    if 'sprot_best' in leaf and not available(leaf['sprot_best']) and best_hit_product(link):
                        raise ValueError('Distribution annotation best hit disagrees with the exact family leaf')
                    predicted = annotation['swissprot_best_hit_protein_name']
                    if predicted and best_hit_product(link) and predicted != best_hit_product(link):
                        raise ValueError('Distribution hit protein name disagrees with the exact family leaf')
                    if predicted:
                        link['swissprot_best_hit_protein_name'] = predicted
                for row in rows:
                    if row.get("child1") == row.get("child2") == "-999":
                        matches = [s for s in species if row["node_name"].startswith(s + "_")]
                        if not matches:
                            raise ValueError("Gene-tree tip cannot be mapped to analysis species: " + row["node_name"])
                        membership[family, max(matches, key=len)] += 1
            for logical, expected in sources.items():
                with store.open_binary(*logical.split("/", 1)) as h:
                    if hashlib.sha256(h.read()).hexdigest() != expected:
                        raise ValueError("Family input changed during distribution plotting")
    annotations.verify()
    from focus_hgt_gene_trees import background_supported

    selected_ids = {r["event_id"] for r in selected}
    labels = product_labels(families, annotation_links)
    highlighted = {
        (r["orthogroup"], r["gene_species"])
        for r in links
        if r["event_id"] in selected_ids and r["side"] == "recipient" and background_supported(r)
    }
    distribution = [
        dict(
            orthogroup=f,
            species=s,
            category1=values.get(s),
            gene_tree_tip_count="" if f in missing else membership[f, s],
            supported_recipient=int((f, s) in highlighted),
            **labels[f],
        )
        for f in families
        for s in species
    ]
    write(
        directory / "orthogroup_species_distribution.tsv",
        distribution,
        [
            "orthogroup",
            "species",
            "category1",
            "gene_tree_tip_count",
            "supported_recipient",
            "protein_product",
            "annotation_gene_id",
            "annotation_basis",
            "all_recipient_product_labels",
        ],
    )
    cmap = ListedColormap(["#f4f5f6", "#bad2e2", "#6a9cbc", BLUE, "#173e5a"])
    cmap.set_bad("#c3c7ca")
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 4.5, 9.5, 100000], 5)
    matrix = np.array([[np.nan if f in missing else membership[f, s] for s in species] for f in families])
    # A separate annotation column leaves both product labels and species cells readable.
    fig = plt.figure(figsize=(23, max(9, len(families) * 0.35 + 4)))
    grid = fig.add_gridspec(
        1, 2, left=0.035, right=0.985, top=0.80, bottom=0.34, width_ratios=[0.33, 0.67], wspace=0.01
    )
    names = fig.add_subplot(grid[0])
    ax = fig.add_subplot(grid[1])
    n = max(1, len(families))
    names.set_xlim(0, 1)
    names.set_ylim(n - 0.5, -0.5)
    names.axis("off")
    names.text(0.01, -1.1, "Orthogroup", fontsize=10, weight="bold")
    names.text(0.19, -1.1, "Protein product / best-hit prediction", fontsize=10, weight="bold")
    for i, f in enumerate(families):
        names.text(0.01, i, f, fontsize=8, va="center", color=BLUE)
        names.text(0.19, i, textwrap.fill(labels[f]["protein_product"], width=53), fontsize=8, va="center")
    if families:
        ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    ax.set_yticks([])
    ax.set_xticks(range(len(species)), [s.replace("_", " ") for s in species], rotation=90, fontsize=7)
    ax.set_ylim(n - 0.5, -0.5)
    ax.tick_params(length=0)
    for j, s in enumerate(species):
        if values.get(s) == 1:
            ax.get_xticklabels()[j].set_color(ORANGE)
            ax.add_patch(Rectangle((j - 0.5, -1.2), 1, 0.4, color=ORANGE, clip_on=False))
    for i, f in enumerate(families):
        for j, s in enumerate(species):
            if (f, s) in highlighted:
                ax.add_patch(Rectangle((j - 0.43, i - 0.43), 0.86, 0.86, fill=False, ec=ORANGE, lw=1.3))
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i, (label, color) in enumerate(zip(["0", "1", "2–4", "5–9", "10+"], cmap.colors, strict=True)):
        fig.patches.append(
            Rectangle((0.55 + i * 0.065, 0.856), 0.013, 0.014, transform=fig.transFigure, facecolor=color)
        )
        fig.text(0.566 + i * 0.065, 0.863, label, fontsize=8, va="center")
    save(
        fig,
        "orthogroup_species_distribution",
        "Where the focused orthogroups occur in analyzed trees",
        f"{len(families)} orthogroups × {len(species)} species in species-tree tip order. Orange outlines: individually supported category-1 recipient genes.",
        "Colors count existing gene-tree tips, including duplicate copies; zero does not establish biological absence. Missing families remain NA.\n"
        "Product labels always use one retained recipient best-hit prediction per family; source gene, evidence basis and all available labels are in the TSV. Predicted labels do not establish function.",
    )
    nodes = {node.name: [t.name for t in node.get_terminals()] for node in tree.find_clades() if node.name}

    def donor_group(event):
        branch = event["generax_donor_node"]
        classes = event.get("donor_descendant_classes", "")
        return (
            classes.replace(";", " + ") + f" ({branch})"
            if branch in nodes and len(nodes[branch]) > 1 and classes
            else classes or branch.replace("_", " ")
        )

    donors = sorted({donor_group(e) for e in selected})
    recipients = sorted({e["generax_recipient_node"] for e in selected})
    pairs = Counter((donor_group(e), e["generax_recipient_node"]) for e in selected)
    mapped = [
        dict(
            event_id=e["event_id"],
            orthogroup=e["orthogroup"],
            donor_group=donor_group(e),
            donor_branch=e["generax_donor_node"],
            recipient_branch=e["generax_recipient_node"],
            donor_clade_tip_labels="; ".join(nodes.get(e["generax_donor_node"], [])),
            recipient_clade_tip_labels="; ".join(nodes.get(e["generax_recipient_node"], [])),
        )
        for e in selected
    ]
    write(
        directory / "donor_recipient_events.tsv",
        mapped,
        [
            "event_id",
            "orthogroup",
            "donor_group",
            "donor_branch",
            "recipient_branch",
            "donor_clade_tip_labels",
            "recipient_clade_tip_labels",
        ],
    )
    write(
        directory / "donor_recipient_counts.tsv",
        [dict(donor_group=d, recipient_branch=r, event_count=pairs[d, r]) for d in donors for r in recipients],
        ["donor_group", "recipient_branch", "event_count"],
    )
    fig, ax = plt.subplots(figsize=(14, 9))
    fig.subplots_adjust(left=0.26, right=0.91, top=0.80, bottom=0.30)
    maximum = max(pairs.values(), default=1)
    if donors and recipients:
        im = ax.imshow(
            [[pairs[d, r] for r in recipients] for d in donors], cmap="Blues", vmin=0, vmax=maximum, aspect="auto"
        )
        for i, d in enumerate(donors):
            for j, r in enumerate(recipients):
                if pairs[d, r]:
                    ax.text(
                        j,
                        i,
                        str(pairs[d, r]),
                        ha="center",
                        va="center",
                        color="white" if pairs[d, r] > maximum * 0.5 else "#183245",
                        fontsize=11,
                    )
        fig.colorbar(im, ax=ax, shrink=0.75, label="Modeled events", ticks=range(maximum + 1), pad=0.02)
    ax.set_yticks(range(len(donors)), donors, fontsize=9)
    recipient_labels = [
        r.replace("_", " ") if len(nodes.get(r, [])) <= 1 else r + ": " + " + ".join(s.split("_")[0] for s in nodes[r])
        for r in recipients
    ]
    ax.set_xticks(range(len(recipients)), recipient_labels, rotation=35, ha="right", fontsize=9)
    ax.tick_params(length=0)
    title = "Donor groups and " + trait + "-recipient branches"
    save(
        fig,
        "donor_recipient_counts",
        title,
        f"{len(selected)} distinct events; {len(families)} orthogroups. Groups come from modeled donor species branches.",
        "Each ancestral recipient event is counted once; descendant species and post-transfer copies do not multiply event counts.\n"
        "Internal branch IDs and complete descendant tip labels are retained in donor_recipient_events.tsv. No ancestral trait state is inferred.",
    )
    return dict(pdf_count=3, filtering_counts=counts, family_source_sha256=sources,
                context_annotation_source_sha256=annotations.sources)

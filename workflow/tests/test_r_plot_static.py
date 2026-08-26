import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SUPPORT_DIR = REPO_ROOT / "workflow" / "support"
NUMERIC_LEGEND_POSITION_RE = re.compile(r"legend\.position\s*=\s*c\(")


def test_r_scripts_do_not_use_deprecated_numeric_legend_position():
    offending = []
    for path in SUPPORT_DIR.rglob("*"):
        if path.suffix not in {".r", ".R"}:
            continue
        text = path.read_text(encoding="utf-8")
        if NUMERIC_LEGEND_POSITION_RE.search(text):
            offending.append(path.relative_to(REPO_ROOT).as_posix())

    assert offending == [], "Found deprecated numeric legend.position usage in: " + ", ".join(offending)


def test_amino_acid_site_panel_noops_on_empty_site_list():
    stat_branch2tree = (SUPPORT_DIR / "stat_branch2tree_plot.r").read_text(encoding="utf-8")
    fimo_motif = (SUPPORT_DIR / "treevis" / "R" / "05_fimo_motif.R").read_text(encoding="utf-8")

    assert "parse_amino_acid_site_list = function(site_text)" in stat_branch2tree
    assert "return(integer(0))" in stat_branch2tree
    assert "if (length(selected_amino_acid_sites) == 0)" in stat_branch2tree
    assert "amino_acid_site column will not be added" in stat_branch2tree
    assert "return(g)" in fimo_motif


def test_query2family_busco_plot_uses_shared_species_label_parser():
    text = (SUPPORT_DIR / "plot_query2family_presence_absence.R").read_text(encoding="utf-8")

    assert 'source(file.path(script_dir, "species_label_utils.r"), local = TRUE)' in text
    assert "gg_species_label_from_filename(x)" in text
    assert 'sub("_sp_[^._-]+$", "_sp", x)' not in text


def test_query2family_ortholog_legend_expands_plot_bounds():
    text = (SUPPORT_DIR / "plot_query2family_presence_absence.R").read_text(encoding="utf-8")

    assert "ortholog_legend_right <- max(" in text
    assert "x_max <- max(x_max, ortholog_legend_right + 0.15)" in text


def test_query2family_plot_supports_query_gene_orthology_glyphs():
    text = (SUPPORT_DIR / "plot_query2family_presence_absence.R").read_text(encoding="utf-8")

    assert 'ortholog_column_path <- get_arg(args, "ortholog_column_table")' in text
    assert 'ortholog_glyph_path <- get_arg(args, "ortholog_glyph_table")' in text
    assert 'ortholog_tree_path <- get_arg(args, "ortholog_tree_table")' in text
    assert 'ortholog_synteny_path <- get_arg(args, "ortholog_synteny_table")' in text
    assert 'ortholog_ufboot_path <- get_arg(args, "ortholog_ufboot_table")' in text
    assert 'species_mapping_tree_path <- get_arg(args, "species_mapping_tree")' in text
    assert 'reference_species <- get_arg(args, "reference_species")' in text
    assert 'evidence_layout <- tolower(get_arg(args, "evidence_layout", "band"))' in text
    assert 'evidence_layout %in% c("band", "rail", "glyph", "off")' in text
    assert 'evidence_strip_mode <- evidence_layout %in% c("band", "rail")' in text
    assert 'glyph_rect_df$relation == "shared_ancestral"' in text
    assert '"#9ecae1"' in text
    assert "label = copy_label" in text
    assert '"anchor_count"' in text
    assert 'viridis_anchor_colors <- function(anchor_count, max_anchor_count)' in text
    assert 'palette = "viridis"' in text
    assert 'viridis_unit_colors <- function(fraction)' in text
    assert 'viridis_anchor_rail_colors <- function(anchor_count, max_anchor_count)' in text
    assert 'viridis_unit_colors(anchor_count / max_anchor_count)' in text
    assert 'colors <- rep("#f7f7f7", length(anchor_count))' not in text
    assert 'FUN = max' in text
    assert 'shape = 23' in text
    assert '"Local synteny anchors"' in text
    assert '"Local synteny anchors (A)"' not in text
    assert '"Color: highest per-copy anchor count"' in text
    assert 'paste0(evidence_layout_label, " (top): local synteny anchors")' in text
    assert '"Fill: highest per-copy anchor count (0 = none)"' in text
    assert '"Color: highest per-copy A in cell"' not in text
    assert '" downstream gene models"' in text
    assert '" gene models per side (family-specific)"' in text
    assert "marker_label" not in text
    assert "query_label_gap <- 0.06" in text
    assert "species_label_matrix_gap <- 0.18" in text
    assert "heatmap_left <- label_right + species_label_matrix_gap" in text
    assert "size = 0.92" in text
    assert 'viridis_ufboot_colors <- function(ufboot)' in text
    assert 'viridis_ufboot_colors(ufboot_evidence_df$ufboot)' in text
    assert 'viridis_ufboot_colors(ufboot_breaks)' in text
    assert 'inferno_ufboot_colors <- function(ufboot)' in text
    assert 'palette = "Inferno"' in text
    assert '"Gene tree UFBoot (%)"' in text
    assert '"GeneRax-topology UFBoot (%)"' not in text
    assert '"ufboot_support_source"' in text
    assert '"Color: lowest per-copy MRCA-branch support"' in text
    assert '"No circle: unavailable or reference self"' in text
    assert 'paste0(evidence_layout_label, " (bottom): ")' in text
    assert '"Numeric fill requires support for every copy"' in text
    assert 'label = paste0(evidence_layout_label, " states")' in text
    assert 'displayed_evidence_statuses <- c("unavailable", "reference_self")' in text
    assert 'length(displayed_evidence_statuses) > 0' in text
    assert '"No band: reference self"' not in text
    assert 'evidence_state_legend_df$draw_swatch <- rep(' in text
    assert 'synteny_legend_visible <- any(' in text
    assert 'ufboot_legend_visible <- any(' in text
    assert (
        'stop("Ortholog synteny reference_self status disagrees with candidate/reference IDs")'
        in text
    )
    assert 'evidence_df$evidence_status != "reference_self"' not in text
    assert "synteny_display_df <- synteny_evidence_df" in text
    assert "ufboot_display_df <- ufboot_evidence_df" in text
    assert "ufboot_evidence_df$evaluated_count == ufboot_evidence_df$pair_count" in text
    assert "ufboot_evidence_df$pair_count != ufboot_evidence_df$glyph_copy_number" in text
    assert 'synteny_evidence_df$pair_count != synteny_evidence_df$glyph_copy_number' in text
    assert 'evidence_strip_mode && nrow(synteny_evidence_df) > 0' in text
    assert 'identical(evidence_layout, "glyph") && nrow(synteny_marker_df) > 0' in text
    assert 'add_evidence_geometry <- function(evidence_df, evidence_half' in text
    assert 'band_height <- (lane_ymax - lane_ymin) * 0.18' in text
    assert 'evidence_df$evidence_xmin <- ifelse(' in text
    assert 'evidence_df$evidence_xmax <- ifelse(' in text
    assert 'glyph_rect_df$x <- (glyph_rect_df$xmin + glyph_rect_df$xmax) / 2' in text
    assert "FUN = min" in text
    assert 'paste0(reference_species_display, " orthologs")' in text
    assert 'paste0(reference_species_display, " genes")' not in text
    assert 'paste0(reference_species_display, " gene tree")' in text
    assert "Query orthology" not in text
    assert '"reference-gene-specific"' in text
    assert '"query-specific"' not in text
    assert 'query_tree_nodes_df$event == "D"' in text
    assert '"D#: mapped duplication"' in text
    assert '"mapped_species_node"' in text
    assert "species_mapping_label_by_node" in text
    assert "transfer_species_mapping_labels <- function(target_tree, mapping_tree)" in text
    assert "ambiguous_mapped_labels" in text
    assert "duplication_map_df" in text
    assert "root_stem_needed" in text
    assert "root_stem_df" in text
    assert "duplication_map_df$branch_start_raw" in text
    assert "minimum_bar_center_spacing <- 0.16" in text
    assert "duplication_map_df$marker_x[collision_rows] <- preferred_centers" in text
    assert "device = grDevices::cairo_pdf" in text
    assert 'label = "Duplication bar colors"' in text
    assert "duplication_family_mode" in text
    assert 'palette = "Dark 3"' in text
    assert "duplication_map_df$bar_width" in text
    assert "duplication_map_df$ymax <- duplication_map_df$marker_y +" in text
    assert "duplication_map_df$duplication_count / duplication_count_scale_max" in text
    assert "nice_count_breaks <- function(max_count, target_intervals = 3)" in text
    assert "breaks <- pretty(c(0, max_count), n = target_intervals)" in text
    assert "duplication_count_scale_max <- max(duplication_count_breaks)" in text
    assert "count_key_values <- duplication_count_breaks[duplication_count_breaks > 0]" in text
    assert "count_key_title_y <- min(duplication_family_key_df$y, na.rm = TRUE) - 0.90" in text
    assert "count_key_y <- count_key_title_y - 0.85" in text
    assert "family_key_two_column_width <= tree_width / 2" in text
    assert '"non-contiguous orthology"' in text
    assert "data = duplication_map_df" in text
    assert "aes(\n        xmin = xmin," in text
    assert "aes(x = x, y = y, fill = family_color)" in text
    assert 'duplication_family_mode & query_tree_nodes_df$event == "D"' in text
    assert "shape = 21" in text
    assert "duplication_track_mode" not in text
    assert "duplication_circle_mode" not in text
    assert "duplication_family_mode <- has_mapped_duplications" in text
    assert "ortholog_family_count > 1" not in text
    assert '"in_reference_tree"' in text
    assert "duplication_count ~ family_id + family_order + mapped_species_node" in text
    assert "Mapped duplication species nodes were not found in the species mapping tree" in text
    assert "scale_size_area(" not in text
    assert 'label = "Bar height = duplication count"' in text
    assert "label = as.character(family_headers$family_id)" in text
    assert 'paste0("F", family_headers$family_order' not in text
    assert 'paste0("F", query_tree_title_df$family_order' not in text
    assert "family_boundary_df" in text
    assert 'query_tree_title_df$family_id' in text
    query_tree_title_layer = text.split("if (nrow(query_tree_title_df) > 0) {", 2)[2].split(
        "combined <- combined +\n  geom_text(data = query_label_df", 1
    )[0]
    assert 'fontface = "plain"' in query_tree_title_layer
    assert 'fontface = "bold"' not in query_tree_title_layer


def test_query2family_tree_axis_uses_bar_by_default():
    text = (SUPPORT_DIR / "plot_query2family_presence_absence.R").read_text(encoding="utf-8")

    assert 'tree_scale_mode <- get_arg(args, "tree_scale", "bar")' in text
    assert 'dated_ruler_mode <- tree_scale_mode == "ruler"' in text
    assert 'stop("--tree_scale=ruler requires an ultrametric species tree")' in text
    assert 'tree_axis_label <- "Branch length"' in text
    assert "branch_scale_width <- nice_scale_width(" in text
    assert "branch_scale_df <- data.frame(" in text
    assert "data = branch_scale_df" in text
    assert "linewidth = 0.75" in text
    assert 'tree_axis_label <- "Million years ago"' in text
    assert "if (dated_ruler_mode)" in text
    assert "label = tree_axis_label" in text


def test_query2family_copy_numbers_are_layered_above_evidence_bands():
    text = (SUPPORT_DIR / "plot_query2family_presence_absence.R").read_text(encoding="utf-8")

    glyph_rect_layer = text.index(
        "geom_rect(data = glyph_rect_df, aes(xmin = xmin, xmax = xmax,"
    )
    synteny_band_layer = text.index("data = synteny_display_df", glyph_rect_layer)
    ufboot_band_layer = text.index("data = ufboot_display_df", synteny_band_layer)
    copy_number_layer = text.index(
        "aes(x = x, y = text_y, label = copy_label, color = text_color)",
        ufboot_band_layer,
    )

    assert glyph_rect_layer < synteny_band_layer < ufboot_band_layer < copy_number_layer


def test_tree_plot_consumers_load_the_installed_treevis_package():
    for name in ("annotation_summary.r", "stat_branch2tree_plot.r"):
        text = (SUPPORT_DIR / name).read_text(encoding="utf-8")
        assert "library(genegalleon.treevis)" in text
        assert "treevis_dir" not in text
        assert "tree_annotation_dir" not in text
        assert "R/main.R" not in text

    assert not (SUPPORT_DIR / "tree_annotation").is_symlink()
    assert not (SUPPORT_DIR / "treevis" / "R" / "main.R").exists()

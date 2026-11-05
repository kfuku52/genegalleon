source("workflow/support/tree_topology_helpers.R")
source("workflow/support/tree_label_formatting.R")
ordered_tree <- ape::read.tree(text="(A,(B,C),(D,E));")
stopifnot(identical(gg_descendant_tip_labels(ordered_tree,6L),c("A","B","C","D","E")))

source_tree <- ape::read.tree(text = "((A:1,B:1)90:1,(C:1,D:1)90:1)100;")
target_tree <- ape::read.tree(text = "((B:1,A:1):1,(D:1,C:1):1);")
result <- gg_transfer_internal_labels(target_tree, source_tree)
stopifnot(any(result$node.label %in% c("80", "90")))
stopifnot(identical(result$tip.label, target_tree$tip.label))
stopifnot(identical(result$edge, target_tree$edge))
stopifnot(identical(result$edge.length, target_tree$edge.length))

subset_tree <- ape::drop.tip(target_tree, "D")
rejected <- suppressWarnings(gg_transfer_internal_labels(subset_tree, source_tree))
stopifnot(identical(rejected, subset_tree))
accepted <- gg_transfer_internal_labels(subset_tree, source_tree, allow_subset = TRUE)
stopifnot(identical(accepted$tip.label, subset_tree$tip.label))
stopifnot(any(accepted$node.label %in% c("80", "90")))
unlabelled <- source_tree
unlabelled$node.label <- NULL
stopifnot(identical(gg_transfer_internal_labels(target_tree, unlabelled), target_tree))

labels <- c("Oryza_sativa_Indica_Group", "Oryza_sativa_tropical_japonica_subgroup",
            "Cenchrus_americanus_x_Cenchrus_purpureus")
for (label in labels) {
  stopifnot(identical(extract_species_label(paste0(label, "_gene1")), label))
  stopifnot(identical(genegalleon.treevis:::treevis_extract_species_label(paste0(label, "_gene1")), label))
}
cat("Shared tree topology and species-label contracts passed.\n")

local({
  previous <- Sys.getenv("TREEVIS_SPECIES_PARSER", unset=NA_character_)
  original_dir <- getwd()
  directory <- tempfile("species-colours-")
  dir.create(directory)
  on.exit({
    setwd(original_dir)
    unlink(directory, recursive=TRUE)
    if (is.na(previous)) Sys.unsetenv("TREEVIS_SPECIES_PARSER")
    else Sys.setenv(TREEVIS_SPECIES_PARSER=previous)
  })
  setwd(directory)
  Sys.setenv(TREEVIS_SPECIES_PARSER="extended-taxonomic")
  plot <- list(data=data.frame(label=paste0(labels,"_gene1"),
                             spnode_coverage=labels, node=seq_along(labels)))
  coloured <- getFromNamespace("append_branch_tiplab_colors", "genegalleon.treevis")(plot,"species",NULL)
  stopifnot(length(unique(coloured$data$tiplab_color)) == length(labels))
  stopifnot(all(coloured$data$tiplab_color != "#000000"))
})

for (text in c("((A:1,B:1)90:1,(C:1,D:1)80:1);", "((C:1,D:1)80:1,(A:1,B:1)90:1);")) {
  failed <- tryCatch({
    gg_transfer_internal_labels(target_tree, ape::read.tree(text = text))
    FALSE
  }, error = function(e) grepl("NWKIT support-label transfer failed", conditionMessage(e), fixed = TRUE))
  stopifnot(failed)
}
deep <- ape::stree(4096, type = "left")
deep$node.label <- rep("90", deep$Nnode)
transferred <- gg_transfer_internal_labels(deep, deep)
stopifnot(identical(transferred$node.label, c("", rep("90", deep$Nnode - 1L))))
stopifnot(identical(transferred$edge, deep$edge), identical(transferred$tip.label, deep$tip.label))

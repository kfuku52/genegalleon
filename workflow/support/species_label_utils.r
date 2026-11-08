# Source the same species-label implementation used by the installed treevis package.
gg_species_helper_sources <- Filter(function(path) is.character(path) && length(path) == 1L && nzchar(path),
                              lapply(sys.frames(), function(frame) frame$ofile))
gg_species_helper_source <- tail(gg_species_helper_sources, 1L)[[1L]]
source(file.path(dirname(normalizePath(gg_species_helper_source, mustWork = TRUE)),
                 "treevis", "R", "00_00_species_labels.R"), local = TRUE)
rm(gg_species_helper_source, gg_species_helper_sources)

# Tree Visualization helper loader
#
# This file keeps compatibility with the legacy entry point while the real
# implementation is split across module files in the same directory.

this_file = sys.frames()[[1]][['ofile']]
if (is.null(this_file) || !nzchar(this_file)) {
    stop('Could not resolve treevis helper directory.')
}
r_dir = dirname(normalizePath(this_file, winslash = '/', mustWork = TRUE))
module_files = c(
    '00_core.R',
    '01_alignment_panels.R',
    '02_domain_ortholog.R',
    '03_ortholog_tree.R',
    '04_branch_synteny.R',
    '05_fimo_motif.R'
)
for (module_file in module_files) {
    source(file.path(r_dir, module_file), local = TRUE)
}

# ggplot2 aesthetics are evaluated in data masks. Declaring their column names
# here lets R CMD check distinguish those columns from unresolved R objects.
utils::globalVariables(c(
    'alpha_value', 'branch', 'branch_id', 'cluster_membership', 'count', 'end',
    'fill', 'group', 'group_id', 'hjust', 'isTip', 'key', 'label', 'left_end',
    'mid_end', 'mid_start', 'motif_altid', 'name', 'nearests', 'node_category',
    'panel_x', 'plot_value', 'polygon_id', 'position', 'right_start', 'sacc',
    'show_circle', 'species', 'start', 'value', 'x', 'x_dummy', 'x_end',
    'x_start', 'xend', 'xmax', 'xmid', 'xmin', 'y', 'y_end', 'y_start', 'yend',
    'ymax', 'ymid', 'ymin'
))

# Suppresses R CMD check's "no visible binding for global variable" NOTE for names that are only ever
# column references inside non-standard-evaluation contexts (aes(), data.frame() built then referenced by
# column name in a later aes() call, etc.) -- they are never meant to be free-standing variables.
utils::globalVariables(c(
  ".data", "auto_cor", "bin", "block", "coord_x", "coord_y", "count", "curve", "density", "gene.list",
  "group", "iteration", "label", "layer", "method", "ndensity", "param", "parameter", "rate",
  "residuals.log", "slope_seg_neg_x", "slope_seg_neg_y", "slope_seg_pos_x", "slope_seg_pos_y", "step",
  "sw_stat", "theoretical", "timeseries", "treatment", "type", "vals", "value", "value_high",
  "value_low", "x", "x0", "x1", "y", "y1", "y_bins"
))

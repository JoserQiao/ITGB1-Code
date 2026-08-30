packages <- c("readr", "dplyr", "tidyr", "ggplot2", "patchwork", "forcats", "scales")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
message("R dependencies are available.")

# Installs the R packages used by R/Fungi_Paper_Script.R.
# Run once before the analysis. The main script does not install anything.

cran_packages <- c(
  "dplyr", "tidyr", "purrr", "tibble", "ggplot2", "rstatix", "ggpubr",
  "patchwork", "Polychrome", "Hmisc", "ggvenn", "circlize", "writexl",
  "vegan", "ggforce", "igraph", "ggraph", "tidygraph", "ggnewscale",
  "scales", "ggrepel", "cowplot", "FSA", "rcompanion",
  # Needed only when run_anova_comparison <- TRUE
  "multcompView", "car",
  # Used to install the GitHub and Bioconductor packages below
  "remotes", "BiocManager"
)
install.packages(setdiff(cran_packages, rownames(installed.packages())))

bioc_packages <- c("phyloseq", "ComplexHeatmap")
BiocManager::install(setdiff(bioc_packages, rownames(installed.packages())))

# Not on CRAN
if (!requireNamespace("SpiecEasi", quietly = TRUE)) {
  remotes::install_github("zdk123/SpiecEasi")
}
if (!requireNamespace("linkET", quietly = TRUE)) {
  remotes::install_github("Hy4m/linkET")
}

# //////////////////////////////////////////////////////////////////////////////

# Soil fungal diversity and association networks in the Bwindi landscape
# Author: Rapheal Wangalwa and collaborators
# Purpose: manuscript analyses, figures and supplementary tables.
# Manuscript: [TITLE TO BE ADDED], submitted to Landscape Ecology.
#
# How to run
#   1. Open R with the repository root as the working directory (for example,
#      by opening the .Rproj file in RStudio).
#   2. Place the input phyloseq object at data/fungi.RData (see data/README.md).
#   3. Run: source("R/Fungi_Paper_Script.R")
#
# Input: an RDS-serialized phyloseq object with an existing Guild_class field.
# The original input is called fungi.RData but was saved for readRDS().
# Paths in the configuration below are relative to the repository root.
# All figures and tables are written to output_dir; filenames are retained.
# Existing files with matching names will be replaced.
#
# Data availability: the input data are not included in this repository and
# will be made publicly available upon publication of the manuscript.
#
# Dependencies must be installed separately; this script does not install them.
# CRAN: dplyr, tidyr, purrr, tibble, ggplot2, rstatix, ggpubr, patchwork,
# Polychrome, Hmisc, ggvenn, circlize, writexl, vegan, ggforce, igraph,
# ggraph, tidygraph, ggnewscale, scales, ggrepel, cowplot, multcompView,
# FSA, rcompanion, car.
# Bioconductor: phyloseq, ComplexHeatmap.
# Other repositories: SpiecEasi; linkET (for the Mantel correlation figure).

# //////////////////////////////////////////////////////////////////////////////

# -----------------------------------------------------------------------------
# 01. Configuration
# -----------------------------------------------------------------------------
input_file <- file.path("data", "fungi.RData")
output_dir <- "outputs"
run_anova_comparison <- FALSE
random_seed <- 123L
clr_pseudocount <- 1
ordination_permutations <- 999L
nmds_trymax <- 100L

run_fungi_analysis <- function() {
  required_packages <- c(
    "phyloseq", "ComplexHeatmap", "Hmisc", "igraph", "SpiecEasi",
    "ggraph", "tidygraph", "ggnewscale", "ggplot2", "rstatix",
    "ggpubr", "patchwork", "Polychrome", "ggvenn", "circlize",
    "writexl", "vegan", "ggforce", "scales", "ggrepel", "cowplot",
    "FSA", "rcompanion", "tibble", "purrr", "tidyr", "dplyr", "grid"
  )
  if (run_anova_comparison) {
    required_packages <- c(required_packages, "multcompView", "car")
  }
  required_packages <- c(required_packages, "linkET")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Install missing packages: ", paste(missing_packages, collapse = ", "))
  }
  for (pkg in required_packages) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  if (!file.exists(input_file)) stop("Input file not found: ", input_file)
  input_path <- normalizePath(input_file, mustWork = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output_dir.")
  previous_dir <- getwd()
  setwd(output_dir)
  on.exit(setwd(previous_dir), add = TRUE)
  set.seed(random_seed)
  fungi_physeq <- readRDS(input_path)
  if (!inherits(fungi_physeq, "phyloseq")) {
    stop("The input must be an RDS-serialized phyloseq object.")
  }
  required_ranks <- c(
    "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species",
    "Guild_class"
  )
  missing_ranks <- setdiff(required_ranks, rank_names(fungi_physeq))
  if (length(missing_ranks)) {
    stop("Missing taxonomy columns: ", paste(missing_ranks, collapse = ", "))
  }
  metadata <- as(sample_data(fungi_physeq), "data.frame")
  if (!"Veg" %in% names(metadata)) stop("Missing habitat column: Veg")
  if (anyNA(metadata$Veg)) stop("Veg contains missing habitat assignments.")
  veg_levels <- if (is.factor(metadata$Veg)) {
    levels(droplevels(metadata$Veg))
  } else {
    c("SMF", "MF", "UMF", "CR")
  }
  if (!setequal(unique(as.character(metadata$Veg)), c("SMF", "MF", "UMF",
        "CR"))) {
    stop("Veg must contain the four habitats: SMF, MF, UMF and CR.")
  }
  sample_data(fungi_physeq)$Veg <- factor(metadata$Veg, levels = veg_levels)
  if (any(sample_sums(fungi_physeq) <= 0)) {
    stop("Zero-read samples must be resolved before analysis.")
  }
  save_pdf <- function(filename, plot, width, height) {
    ggsave(filename, plot, width = width, height = height,
      bg = "white", device = cairo_pdf)
  }

  # ----------------------------------------------------------------------------
  # 02. Sequencing overview
  # ----------------------------------------------------------------------------
  # Sequencing read summary — standard "reads overview" reported at the start
  # of Results, before diversity/composition analyses

  # 1. Basic dataset dimensions
  n_samples <- nsamples(fungi_physeq)
  n_taxa <- ntaxa(fungi_physeq)

  cat("=== Dataset overview ===\n")
  cat("Number of samples:", n_samples, "\n")
  cat("Number of OTUs (fungal taxa):", n_taxa, "\n\n")

  # 2. Reads per sample — total, min, max, mean, median, SD
  read_counts <- sample_sums(fungi_physeq)

  read_summary <- tibble(
    Total_reads = sum(read_counts),
    Min_reads_per_sample = min(read_counts),
    Max_reads_per_sample = max(read_counts),
    Mean_reads_per_sample = round(mean(read_counts), 1),
    Median_reads_per_sample = median(read_counts),
    SD_reads_per_sample = round(sd(read_counts), 1)
  )

  cat("=== Sequencing read summary (final, decontaminated dataset) ===\n")
  print(read_summary)

  # 3. Per-sample read counts, in full — useful supplementary table
  per_sample_reads <- tibble(
    Sample = sample_names(fungi_physeq),
    Veg = as(sample_data(fungi_physeq), "data.frame")$Veg,
    Reads = read_counts
  ) %>%
    dplyr::arrange(desc(Reads))

  cat("\n=== Per-sample read counts ===\n")
  print(per_sample_reads, n = Inf)

  write.csv(per_sample_reads, "SuppTable_reads_per_sample.csv", row.names =
      FALSE)

  # 4. Reads and taxa summarized per habitat type
  reads_by_habitat <- per_sample_reads %>%
    dplyr::group_by(Veg) %>%
    dplyr::summarise(
      N_samples = n(),
      Total_reads = sum(Reads),
      Min_reads = min(Reads),
      Max_reads = max(Reads),
      Mean_reads = round(mean(Reads), 1),
      Median_reads = median(Reads)
    )

  cat("\n=== Reads summarized by habitat type ===\n")
  print(reads_by_habitat)
  write.csv(reads_by_habitat, "Table_reads_by_habitat.csv", row.names = FALSE)

  # 5. OTU/taxon richness per sample (before any diversity analysis) —
  # also commonly reported alongside read counts
  otu_mat <- as(otu_table(fungi_physeq), "matrix")
  if (!taxa_are_rows(fungi_physeq)) otu_mat <- t(otu_mat)

  otus_per_sample <- colSums(otu_mat > 0)

  cat("\n=== OTUs detected per sample ===\n")
  cat("Min:", min(otus_per_sample), "| Max:", max(otus_per_sample),
    "| Mean:", round(mean(otus_per_sample), 1), "| Median:", median(
      otus_per_sample), "\n")

  # 6. Taxonomic assignment summary — how many OTUs resolved at each rank
  # (you already have this from your earlier taxonomic resolution table —
  # reproduced here in condensed form for the reads-overview paragraph)
  tax_df <- as.data.frame(tax_table(fungi_physeq))
  resolved_summary <- sapply(c("Phylum","Class","Order","Family","Genus",
      "Species"), function(rank) {
      sum(!is.na(tax_df[[rank]]) & tax_df[[rank]] != "" & tax_df[[rank]] !=
        "Unclassified")
  })

  cat("\n=== OTUs resolved at each taxonomic rank (whole dataset) ===\n")
  print(resolved_summary)

  # ----------------------------------------------------------------------------
  # 03. Sampling coordinates
  # ----------------------------------------------------------------------------
  # Extract sampling site coordinates directly from fungi_physeq

  # 1. Check what columns are actually available in sample_data first —
  # confirms the exact column names before trying to select them
  cat("Columns available in sample_data(fungi_physeq):\n")
  print(colnames(sample_data(fungi_physeq)))

  # 2. Extract Sample, Veg, Latitude, Longitude (adjust column names if the
  # print above shows something different, e.g. "Lat"/"Long" instead)
  coords_df <- as(sample_data(fungi_physeq), "data.frame") %>%
    rownames_to_column("Sample") %>%
    dplyr::select(Sample, Veg, Latitude, Longitude)
  if (!all(vapply(coords_df[c("Latitude", "Longitude")], is.numeric,
        logical(1)))) {
    stop("Latitude and Longitude must be numeric.")
  }
  cat("\n=== Sampling site coordinates ===\n")
  print(as.data.frame(coords_df), row.names = FALSE)
  write.csv(coords_df, "sampling_site_coordinates.csv", row.names = FALSE)

  # 3. Quick summary — coordinate range per habitat type, useful for
  # sanity-checking and for setting map extent later
  coords_summary <- coords_df %>%
    dplyr::group_by(Veg) %>%
    dplyr::summarise(
      n = n(),
      Lat_min = min(Latitude, na.rm = TRUE),
      Lat_max = max(Latitude, na.rm = TRUE),
      Lon_min = min(Longitude, na.rm = TRUE),
      Lon_max = max(Longitude, na.rm = TRUE)
    )

  cat("\n=== Coordinate range summary by habitat type ===\n")
  print(coords_summary)

  # ----------------------------------------------------------------------------
  # 04. Environmental contrasts
  # ----------------------------------------------------------------------------
  # 1. Environmental data — numeric variables, same set used throughout
  # (adjust column selection if you want a specific subset)
  env_df_full <- as(sample_data(fungi_physeq), "data.frame") %>%
    rownames_to_column("Sample") %>%
    dplyr::select(Sample, Veg, where(is.numeric), -any_of(c("Latitude",
          "Longitude", "MO")))

  env_vars <- setdiff(colnames(env_df_full), c("Sample", "Veg"))

  if (run_anova_comparison) {
    # 2. Function: for one environmental variable, compute mean±SD per group,
    # run one-way ANOVA, Tukey HSD, and derive compact letter display
    analyze_env_variable <- function(var_name, data, veg_levels) {

      data_sub <- data %>% dplyr::select(Veg, value = all_of(var_name)) %>%
        dplyr::filter(!is.na(value))

      # --- Descriptive stats ---
      desc_stats <- data_sub %>%
        dplyr::group_by(Veg) %>%
        dplyr::summarise(Mean = mean(value), SD = sd(value), N = n(), .groups =
            "drop")

      # --- One-way ANOVA ---
      model <- aov(value ~ Veg, data = data_sub)
      anova_summary <- summary(model)[[1]]
      f_stat <- anova_summary$`F value`[1]
      p_value <- anova_summary$`Pr(>F)`[1]

      # --- Tukey HSD post-hoc + compact letter display ---
      # (only meaningful/stable if the omnibus ANOVA is run; letters are
      # still generated regardless, consistent with standard convention of
      # reporting them alongside the omnibus test)
      tukey <- TukeyHSD(model)
      tukey_pvals <- tukey$Veg[, "p adj"]
      letters_df <- multcompLetters(tukey_pvals)$Letters
      letters_tbl <- tibble(Veg = names(letters_df), Letter = letters_df)

      # --- Combine ---
      result <- desc_stats %>%
        dplyr::left_join(letters_tbl, by = "Veg") %>%
        dplyr::mutate(
          Variable = var_name,
          F_stat = round(f_stat, 3),
          p_value = round(p_value, 4),
          Display = paste0(round(Mean, 2), " \u00b1 ", round(SD, 2),
            letters_tbl$Letter[match(Veg, letters_tbl$Veg)])
        )

      result
    }

    # 3. Run across all environmental variables
    all_results <- map_dfr(env_vars, function(v) {
        tryCatch(
          analyze_env_variable(v, env_df_full, veg_levels),
          error = function(e) {
            cat("Skipped", v, "- error:", conditionMessage(e), "\n")
            NULL
          }
        )
    })

    # 4. Reshape into the final wide table: one row per variable, one column
    # per Veg group showing "Mean ± SD^letter", plus F and p-value columns
    final_table <- all_results %>%
      dplyr::select(Variable, Veg, Display, F_stat, p_value) %>%
      pivot_wider(names_from = Veg, values_from = Display) %>%
      dplyr::distinct(Variable, .keep_all = TRUE) %>%
      dplyr::select(Variable, all_of(veg_levels), F_stat, p_value) %>%
      dplyr::arrange(p_value)

    cat(paste0("=== Descriptive statistics + ANOVA + ",
        "Tukey HSD (compact letters) ===\n"))
    print(final_table, n = Inf)

    write.csv(final_table, "Table_env_variables_ANOVA_TukeyHSD.csv", row.names =
        FALSE)

    # 5. Also save the FULL long-format results (raw group means/SD/N +
    # F/p, without the combined display string) — useful for double-checking
    # or building your own custom table layout later
    write.csv(all_results, "Table_env_variables_full_results_long.csv",
      row.names = FALSE)

    # Optional: flag variables where ANOVA assumptions may be questionable
    check_assumptions <- map_dfr(env_vars, function(v) {
        data_sub <- env_df_full %>% dplyr::select(Veg, value = all_of(v)) %>%
          dplyr::filter(!is.na(value))
        model <- aov(value ~ Veg, data = data_sub)
        shapiro_p <- shapiro.test(residuals(model))$p.value
        levene_p <- car::leveneTest(value ~ Veg, data = data_sub)$`Pr(>F)`[1]
        tibble(Variable = v, Shapiro_p = round(shapiro_p, 4), Levene_p = round(
            levene_p, 4))
    })

    cat(paste0("\n=== Variables with potential ",
        "normality/homogeneity violations (p<0.05) ", "===\n"))
    print(check_assumptions %>% dplyr::filter(Shapiro_p < 0.05 | Levene_p <
        0.05))

  }
  # post-hoc with compact letter display — consistent with the Hill-number
  # analysis, given majority of environmental variables violate ANOVA
  # assumptions (Shapiro/Levene tests)

  analyze_env_variable_np <- function(var_name, data) {

    data_sub <- data %>% dplyr::select(Veg, value = all_of(var_name)) %>%
      dplyr::filter(!is.na(value))

    # --- Descriptive stats: Median and IQR, consistent with non-parametric test
    desc_stats <- data_sub %>%
      dplyr::group_by(Veg) %>%
      dplyr::summarise(
        Median = median(value),
        Q1 = quantile(value, 0.25),
        Q3 = quantile(value, 0.75),
        N = n(),
        .groups = "drop"
      )

    kw <- kruskal.test(value ~ Veg, data = data_sub)

    dunn <- dunnTest(value ~ Veg, data = data_sub, method = "bh")
    dunn_res <- dunn$res

    letters_df <- rcompanion::cldList(P.adj ~ Comparison, data = dunn_res,
      threshold = 0.05)
    letters_tbl <- tibble(Veg = letters_df$Group, Letter = letters_df$Letter)

    desc_stats %>%
      dplyr::left_join(letters_tbl, by = "Veg") %>%
      dplyr::mutate(
        Variable = var_name,
        KW_stat = round(kw$statistic, 3),
        p_value = round(kw$p.value, 4),
        Display = paste0(round(Median, 2), " (", round(Q1, 2), "\u2013", round(
            Q3, 2), ")", Letter)
      )
  }

  all_results_np <- map_dfr(env_vars, function(v) {
      tryCatch(
        analyze_env_variable_np(v, env_df_full),
        error = function(e) {
          cat("Skipped", v, "- error:", conditionMessage(e), "\n")
          NULL
        }
      )
  })

  final_table_np <- all_results_np %>%
    dplyr::select(Variable, Veg, Display, KW_stat, p_value) %>%
    pivot_wider(names_from = Veg, values_from = Display) %>%
    dplyr::distinct(Variable, .keep_all = TRUE) %>%
    dplyr::select(Variable, all_of(veg_levels), KW_stat, p_value) %>%
    dplyr::arrange(p_value)

  cat(paste0("=== Descriptive statistics (Median [IQR]) ",
      "+ Kruskal-Wallis + Dunn's test ===\n"))
  print(final_table_np, n = Inf)

  write.csv(final_table_np,
    "Table_env_variables_KruskalWallis_Dunn_MedianIQR.csv",
    row.names = FALSE, fileEncoding = "UTF-8")

  # ----------------------------------------------------------------------------
  # 05. Taxonomic summaries
  # ----------------------------------------------------------------------------
  # ------------------------------ DOMAIN SUMMARY ------------------------------
  tax_guilded <- as.data.frame(tax_table(fungi_physeq)) %>%
    rownames_to_column("OTU_ID")

  domain_summary <- tax_guilded %>%
    dplyr::count(Domain, name = "OTUs") %>%
    dplyr::mutate(
      Rank = "Domain",
      Percent = 100 * OTUs / sum(OTUs)
    ) %>%
    dplyr::rename(Category = Domain) %>%
    dplyr::select(Rank, Category, OTUs, Percent)

  # -------------------------------- FUNGI SUBSET ------------------------------

  fungi_tax <- tax_guilded %>%
    dplyr::filter(Domain == "Fungi")

  # ---------------------------- OTHER TAXONOMIC SUMMARIES----------------------

  tax_columns <- c(
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "Guild_class"
  )

  summary_list <- lapply(tax_columns, function(col) {

      fungi_tax %>%
        dplyr::count(!!sym(col), name = "OTUs") %>%
        dplyr::mutate(
          Rank = col,
          Percent = 100 * OTUs / sum(OTUs)
        ) %>%
        dplyr::rename(Category = !!sym(col)) %>%
        dplyr::select(Rank, Category, OTUs, Percent)

  })

  fungi_tax_summary <- dplyr::bind_rows(summary_list)

  # -------------------------------- MERGE EVERYTHING --------------------------

  final_summary <- dplyr::bind_rows(
    domain_summary,
    fungi_tax_summary
  )

  # Arrange nicely
  final_summary <- final_summary %>%
    dplyr::arrange(Rank, desc(OTUs))

  # print(final_summary)
  # export the final summary
  write.csv(
    final_summary,
    "Domain_Fungi_OTU_Proportions_New.csv",
    row.names = FALSE
  )
  # Rearrange the table invert wide format with each rank as a separate block of
  # columns
  
  # -------------------------------- FUNCTION ----------------------------------

  make_rank_block <- function(data, rank_name, display_name) {

    data %>%
      dplyr::filter(Rank == rank_name) %>%
      dplyr::arrange(desc(OTUs)) %>%
      dplyr::mutate(
        Percent = round(Percent, 2)
      ) %>%
      dplyr::select(Category, OTUs, Percent) %>%
      dplyr::rename(
        !!display_name := Category,
        !!paste0(display_name, "_OTUs") := OTUs,
        !!paste0(display_name, "_Percent") := Percent
      )

  }

  # ----------------------------------- CREATE BLOCKS --------------------------

  domain_block  <- make_rank_block(final_summary, "Domain", "Domain")

  phylum_block  <- make_rank_block(final_summary, "Phylum", "Phylum")

  class_block   <- make_rank_block(final_summary, "Class", "Class")

  order_block   <- make_rank_block(final_summary, "Order", "Order")

  family_block   <- make_rank_block(final_summary, "Family", "Family")

  genus_block   <- make_rank_block(final_summary, "Genus", "Genus")

  species_block <- make_rank_block(final_summary, "Species", "Species")

  guild_block   <- make_rank_block(final_summary,
    "Guild_class",
    "Trophic_Guild")

  # pad the rows of each block to ensure they all have the same number of rows

  pad_rows <- function(df, n) {

    if (nrow(df) < n) {
      df[(nrow(df) + 1):n, ] <- NA
    }

    df
  }

  max_rows <- max(
    nrow(domain_block),
    nrow(phylum_block),
    nrow(class_block),
    nrow(order_block),
    nrow(family_block),
    nrow(genus_block),
    nrow(species_block),
    nrow(guild_block)
  )
  # make the final table with all the blocks side by side
  final_summary_wide <- dplyr::bind_cols(

    pad_rows(domain_block,  max_rows),

    pad_rows(phylum_block,  max_rows),

    pad_rows(class_block,   max_rows),

    pad_rows(order_block,   max_rows),

    pad_rows(family_block,   max_rows),

    pad_rows(genus_block,   max_rows),

    pad_rows(species_block, max_rows),

    pad_rows(guild_block,   max_rows)

  )
  # print final_summary_wide tbale
  final_summary_wide
  # export the final summary
  write_xlsx(
    final_summary_wide,
    "Final_Summary_Wide_Table_New.xlsx"
  )

  # Taxon-Environment Correlation Heatmap (family-level, CLR-transformed)

  # ----------------------------------------------------------------------------
  # 06. Family and environment associations
  # ----------------------------------------------------------------------------
  set.seed(123)

  # 1. Aggregate to family level

  fungi_family <- tax_glom(
    fungi_physeq,
    taxrank = "Family",
    NArm = FALSE
  )

  family_labels <- as.data.frame(tax_table(fungi_family))$Family

  taxa_names(fungi_family) <- make.unique(as.character(family_labels))

  # Drop unresolved families
  keep_family <- !is.na(taxa_names(fungi_family)) &
    taxa_names(fungi_family) != "" &
    !grepl("^\\?", taxa_names(fungi_family)) &
    !grepl("Unclassified", taxa_names(fungi_family), ignore.case = TRUE)

  fungi_family <- prune_taxa(keep_family, fungi_family)

  cat("Families retained for heatmap:", ntaxa(fungi_family), "\n")

  # 2. CLR-transform family-level abundances

  otu_mat <- as(otu_table(fungi_family), "matrix")

  if (taxa_are_rows(fungi_family)) {
    otu_mat <- t(otu_mat)
  }

  otu_pseudo <- otu_mat + 1

  clr_mat <- t(
    apply(
      otu_pseudo,
      1,
      function(x) log(x) - mean(log(x))
    )
  )

  rownames(clr_mat) <- rownames(otu_mat)

  # 3. Pull environmental variables and align samples

  env_df <- as(sample_data(fungi_family), "data.frame")

  numeric_vars <- env_df %>%
    dplyr::select(where(is.numeric)) %>%
    dplyr::select(-any_of(c("Latitude", "Longitude", "MO")))

  var_check <- sapply(
    numeric_vars,
    function(x) var(x, na.rm = TRUE)
  )

  numeric_vars <- numeric_vars[, var_check > 1e-6, drop = FALSE]

  common_samples <- intersect(
    rownames(clr_mat),
    rownames(numeric_vars)
  )

  clr_mat <- clr_mat[common_samples, , drop = FALSE]
  numeric_vars <- numeric_vars[common_samples, , drop = FALSE]

  # 4. Spearman correlation and FDR correction

  combined_mat <- as.matrix(
    cbind(clr_mat, numeric_vars)
  )

  corr_res <- Hmisc::rcorr(
    combined_mat,
    type = "spearman"
  )

  family_ids <- colnames(clr_mat)
  env_ids <- colnames(numeric_vars)

  r_mat <- corr_res$r[family_ids, env_ids, drop = FALSE]
  p_mat <- corr_res$P[family_ids, env_ids, drop = FALSE]

  p_adj_mat <- matrix(
    p.adjust(as.vector(p_mat), method = "BH"),
    nrow = nrow(p_mat),
    ncol = ncol(p_mat),
    dimnames = dimnames(p_mat)
  )

  # 5. Keep families with at least one significant association

  sig_rows <- rownames(p_adj_mat)[
    apply(
      p_adj_mat,
      1,
      function(x) any(x < 0.05, na.rm = TRUE)
    )
  ]
  # # sig_rows
  # abs(r_mat[sig_rows, , drop = FALSE]),
  # 1,
  # max,
  # )
  #
  # sort(max_r, decreasing = TRUE)
  # )[1:178]
  # view number of families that met the significance threshold
  length(sig_rows)

  r_mat_sig <- r_mat[sig_rows, , drop = FALSE]
  p_adj_sig <- p_adj_mat[sig_rows, , drop = FALSE]

  cat("Families shown in heatmap:", nrow(r_mat_sig), "\n")

  # 6. Build guild annotation for family rows

  guild_lookup <- tax_guilded %>%
    dplyr::filter(!is.na(Family)) %>%
    dplyr::group_by(Family) %>%
    dplyr::summarise(
      Guild_class = names(sort(table(Guild_class), decreasing = TRUE))[1],
      .groups = "drop"
    ) %>%
    dplyr::filter(Family %in% sig_rows)

  row_guild <- guild_lookup$Guild_class[
    match(rownames(r_mat_sig), guild_lookup$Family)
  ]

  row_guild[is.na(row_guild)] <- "Unassigned"
  guild_levels <- unique(row_guild)

  guild_levels
  guild_colors <- c(
    # brown
    "Dung saprotroph"       = "#8B4513",
    # red
    "Fungal parasite"       = "yellow",
    # green
    "Generalist saprotroph" = "darkgreen",
    # purple
    "Lichenized"            = "#FB53ED",
    # orange
    "Mixed"                 = "darkblue",
    # blue
    "Mycorrhizal"           = "#377EB8",
    # dark orange
    "Plant pathogen"        = "#E41A1C",
    # olive green
    "Plant saprotroph"      = "#09BC15",
    # grey
    "Unassigned"               = "#BDBDBD",
    # dark brown
    "Wood saprotroph"       = "#7DFBF0"
  )

  # Ensure colours match only guilds present in your data
  guild_colors <- guild_colors[guild_levels]

  row_anno <- rowAnnotation(
    Guild = row_guild,
    col = list(Guild = guild_colors),
    annotation_legend_param = list(
      Guild = list(title = "Dominant functional guild")
    )
  )

  # 7. Significance stars

  sig_stars <- matrix(
    "",
    nrow = nrow(p_adj_sig),
    ncol = ncol(p_adj_sig),
    dimnames = dimnames(p_adj_sig)
  )

  sig_stars[p_adj_sig < 0.001] <- "***"
  sig_stars[p_adj_sig >= 0.001 & p_adj_sig < 0.01] <- "**"
  sig_stars[p_adj_sig >= 0.01 & p_adj_sig < 0.05] <- "*"

  # 8. Build heatmap

  col_fun <- colorRamp2(
    c(-1, 0, 1),
    c("blue", "white", "red")
  )

  ht <- Heatmap(
    r_mat_sig,
    name = "Spearman r",
    col = col_fun,
    left_annotation = row_anno,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 9),

    layer_fun = function(j, i, x, y, width, height, fill) {
      grid.text(
        sig_stars[cbind(i, j)],
        x,
        y,
        gp = gpar(fontsize = 9, col = "black")
      )
    },

    heatmap_legend_param = list(
      title = "Spearman\nr",
      at = c(-1, -0.5, 0, 0.5, 1)
    ),
    row_title = "Fungal family",
    column_title = "Environmental variable",
    border = TRUE
  )
  ht_drawn <- draw(ht)
  # 9. Export figure and tables
  # Export image
  cairo_pdf(
    "family_environment_heatmap4.pdf",
    width = 17,
    height = 22
  )

  draw(ht_drawn)
  dev.off()

  png(
    "family_environment_heatmap4.png",
    width = 17,
    height = 22,
    units = "in",
    res = 600
  )

  draw(ht_drawn)
  dev.off()
  # Export tables
  r_df <- as.data.frame(r_mat_sig) %>%
    rownames_to_column("Family")

  p_df <- as.data.frame(p_adj_sig) %>%
    rownames_to_column("Family")

  write.csv(
    r_df,
    "family_environment_correlations_r.csv",
    row.names = FALSE
  )

  write.csv(
    p_df,
    "family_environment_correlations_padj.csv",
    row.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # 07. Taxonomic resolution and community overview
  # ----------------------------------------------------------------------------
  # Results table: Resolved taxa, Resolved OTUs, Unresolved OTUs, Total OTUs
  # by Rank and Vegetation type (Veg)

  ranks_to_summarize <- c("Phylum", "Class", "Order", "Family", "Genus",
    "Species")
  veg_groups <- levels(sample_data(fungi_physeq)$Veg)

  # 2. Build the results table
  summarize_group <- function(veg_level) {
    samples_in_group <- sample_names(fungi_physeq)[sample_data(
        fungi_physeq)$Veg == veg_level]
    ps_sub <- prune_samples(samples_in_group, fungi_physeq)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)

    tax_df <- as.data.frame(tax_table(ps_sub))
    total_otus <- ntaxa(ps_sub)

    map_dfr(ranks_to_summarize, function(rnk) {
        vals <- as.character(tax_df[[rnk]])

        is_unresolved <- is.na(vals) | vals == "Unclassified"
        named_vals    <- vals[!is_unresolved]

        # distinct named taxa
        n_resolved_taxa <- length(unique(named_vals))
        # OTUs assigned to a named taxon
        n_resolved_otus <- sum(!is_unresolved)
        # OTUs with no resolution at this rank
        n_unresolved_otus <- sum(is_unresolved)

        tibble(
          Veg = veg_level,
          Rank = rnk,
          `No. of Resolved Taxa` = n_resolved_taxa,
          `Resolved OTUs` = n_resolved_otus,
          `Unresolved (OTUs)` = n_unresolved_otus,
          Total_OTUs = total_otus
        )
    })
  }

  results_table <- map_dfr(veg_groups, summarize_group) %>%
    dplyr::mutate(Rank = factor(Rank, levels = ranks_to_summarize),
      Veg  = factor(Veg, levels = veg_groups)) %>%
    dplyr::arrange(Veg, Rank)

  # 3. Sanity check — Resolved OTUs + Unresolved OTUs must equal Total_OTUs
  stopifnot(all(results_table$`Resolved OTUs` +
        results_table$`Unresolved (OTUs)` == results_table$Total_OTUs))

  cat("=== Results table ===\n")
  print(results_table, n = Inf)

  write.csv(results_table, "taxonomic_resolution_results_table.csv", row.names =
      FALSE)
  # Single combined figure: OTU resolution (stacked, honest scale) +
  # distinct resolved taxa count (text annotation) across Rank, faceted by Veg
  # 1. Reshape results_table to long format for stacking Resolved/Unresolved
  plot_df <- results_table %>%
    pivot_longer(
      cols = c(`Resolved OTUs`, `Unresolved (OTUs)`),
      names_to = "Status", values_to = "N_OTUs"
    ) %>%
    dplyr::mutate(
      Status = factor(Status, levels = c("Unresolved (OTUs)", "Resolved OTUs"),
        labels = c("Unresolved", "Resolved")),
      Rank = factor(Rank, levels = c("Phylum", "Class", "Order", "Family",
          "Genus", "Species")),
      Veg = factor(Veg, levels = levels(results_table$Veg))
    )

  # Label positions: place the distinct-taxa count just above each full bar
  label_df <- results_table %>%
    dplyr::mutate(Rank = factor(Rank, levels = c("Phylum", "Class", "Order",
          "Family", "Genus", "Species")),
      Veg = factor(Veg, levels = levels(results_table$Veg)))

  # 2. Build the combined figure
  status_colors <- c(Resolved = "#6E0065", Unresolved = "#F2B705")

  p_combined <- ggplot(plot_df, aes(x = Rank, y = N_OTUs, fill = Status)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.2) +
    geom_text(
      data = label_df,
      aes(x = Rank, y = Total_OTUs, label = paste0("(", `No. of Resolved Taxa`,
          ")")),
      inherit.aes = FALSE, vjust = -1.8, size = 2, fontface = "bold", color =
        "black"
    ) +
    facet_wrap(~ Veg, nrow = 1) +
    scale_fill_manual(values = status_colors) +
  # more headroom for the label
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(
      x = "Taxonomic Rank", y = "Number of OTUs", fill = NULL,
      title = "Taxonomic resolution across ranks and vegetation types",
      subtitle = paste0("Bars: resolved vs. unresolved OTUs ",
        "(stacked to total). Labels above bars ( ",
        "): distinct resolved taxa count at that ", "rank.")
    ) +
    theme_minimal(base_size = 12.5) +
    theme(
      strip.text = element_text(face = "bold", size = 13),
      strip.background = element_rect(fill = "grey92", color = "grey30",
        linewidth = 0.6),
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.7),
      panel.spacing = unit(0.9, "lines"),
      axis.title.x = element_text(face = "bold", size = 12, margin = margin(t =
            10)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9.5, face =
          "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9.5, color = "black")
    )

  p_combined
  ggsave("taxonomic_resolution_combined.pdf", p_combined, width = 13, height =
      6)
  ggsave("taxonomic_resolution_combined.png", p_combined, width = 13, height =
      6, dpi = 300)

  # Functional guild relative abundance stacked bar chart across Veg
  # 1. Convert to relative abundance BEFORE aggregating by guild
  fungi_rel <- transform_sample_counts(fungi_physeq, function(x) x / sum(x))

  # 2. Melt to long format and sum OTU-level relative abundance within each
  # Guild_class, per sample (an OTU's guild assignment comes from FUNGuild,
  # already stored in tax_table as "Guild_class")
  guild_df <- psmelt(fungi_rel) %>%
    dplyr::mutate(Guild_class = ifelse(is.na(Guild_class) | Guild_class == "",
        "Unknown", Guild_class)) %>%
    dplyr::group_by(Sample, Veg, Guild_class) %>%
    dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")

  # 3. Average relative abundance per Guild_class within each Veg group
  guild_summary <- guild_df %>%
    dplyr::group_by(Veg, Guild_class) %>%
    dplyr::summarise(mean_RA = mean(Abundance), .groups = "drop")

  # Order guilds by overall mean abundance (largest first), "Unknown" always
  # last
  guild_order <- guild_summary %>%
    dplyr::group_by(Guild_class) %>%
    dplyr::summarise(overall_mean = mean(mean_RA)) %>%
    dplyr::arrange(desc(overall_mean)) %>%
    dplyr::pull(Guild_class)

  guild_order <- c(setdiff(guild_order, "Unknown"), "Unknown")

  guild_summary <- guild_summary %>%
    dplyr::mutate(Guild_class = factor(Guild_class, levels = guild_order),
      Veg = factor(Veg, levels = levels(sample_data(fungi_physeq)$Veg)))
  # save guild_summary csv
  write.csv(guild_summary, "guild_summary.csv", row.names = FALSE)

  # 4. Color palette — distinguishable colors for each guild, "Unknown" fixed
  # grey
  # Manual color palette for functional guilds — grouped by color family,
  # spaced apart so adjacent legend entries don't look similar
  guild_colors <- c(
    # blue
    "Mixed"                  = "#0B0A4A",
    # orange
    "Generalist saprotroph"  = "#EF8A17",
    # green
    "Endophyte"               = "green",
    # red
    "Plant pathogen"          = "blue",
    # purple
    "Wood saprotroph"         = "#4AFFF0",
    # brown
    "Dung saprotroph"         = "#8C6D31",
    # mustard/gold
    "Plant saprotroph"        = "yellow",
    # teal-green
    "Mycorrhizal"             = "darkgreen",
    # magenta/pink
    "Fungal parasite"         = "pink",
    # dusty rose
    "Animal pathogen"         = "red",
    # cyan
    "Lichen parasite"         = "#17BECF",
    # yellow-green
    "Lichenized"              = "#B5CF6B",
    # deep plum
    "Nematophagous"           = "#8B38FF",
    # dark brick red
    "Generalist pathogen"     = "#760914",
    # fixed neutral, unassigned
    "Unassigned"              = "grey60"
  )

  p_guild <- ggplot(guild_summary, aes(x = Veg, y = mean_RA, fill =
        Guild_class)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = guild_colors) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(
      x = "Vegetation type", y = "Mean relative abundance", fill =
        "Functional guild",
      title = "Functional guild composition across vegetation types",
      subtitle = paste0("Mean relative abundance per guild ",
        "(FUNGuild assignment), averaged across ", "samples within each group")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(face = "bold", size = 12),
      axis.title.x = element_text(face = "bold", size = 12, margin = margin(t =
            10)),
      legend.position = "right",
      legend.text = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9.5, color = "grey30")
    )

  p_guild
  ggsave("guild_composition_by_Veg.pdf", p_guild, width = 10, height = 7)
  ggsave("guild_composition_by_Veg.png", p_guild, width = 10, height = 7, dpi =
      300)

  # 6. Export the underlying summary table too — useful as a supplementary table
  write.csv(guild_summary, "guild_composition_summary.csv", row.names = FALSE)

  # Venn diagram of shared/unique OTUs across Vegetation types (Veg)
  # 1. Build a named list of OTU IDs present (nonzero abundance) in each
  # Veg group — pooled across all samples within that group
  veg_groups <- levels(sample_data(fungi_physeq)$Veg)

  otu_sets <- lapply(veg_groups, function(v) {
      samples_in_group <- sample_names(fungi_physeq)[sample_data(
          fungi_physeq)$Veg == v]
      ps_sub <- prune_samples(samples_in_group, fungi_physeq)
      ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
      taxa_names(ps_sub)
  })
  names(otu_sets) <- veg_groups

  cat("OTUs detected per group:\n")
  print(sapply(otu_sets, length))

  # 2. Plot the 4-way Venn diagram
  # SMF, MF, UMF, CR
  venn_colors <- c("red", "blue", "green", "orange")

  p_venn <- ggvenn(
    otu_sets,
    fill_color = venn_colors,
    fill_alpha = 0.6,
    stroke_color = "black",
    stroke_size = 0.9,
    set_name_size = 7,
    text_size = 4.5
  ) +
    labs(title = "Shared and unique OTUs across vegetation types") +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

  p_venn

  ggsave("otu_venn_by_Veg.pdf", p_venn, width = 8, height = 8)
  ggsave("otu_venn_by_Veg.png", p_venn, width = 8, height = 8, dpi = 300)

  # 3. Export the actual OTU membership per intersection — useful as a
  # supplementary table (e.g., "OTUs unique to Cropland")
  all_otus <- unique(unlist(otu_sets))

  membership_df <- data.frame(
    OTU = all_otus,
    SMF = all_otus %in% otu_sets[["SMF"]],
    MF  = all_otus %in% otu_sets[["MF"]],
    UMF = all_otus %in% otu_sets[["UMF"]],
    CR  = all_otus %in% otu_sets[["CR"]]
  )

  write.csv(membership_df, "otu_venn_membership.csv", row.names = FALSE)

  # Quick summary: OTUs unique to each group, and core OTUs shared by all 4
  cat("\nOTUs unique to each group:\n")
  print(colSums(membership_df[,-1] & rowSums(membership_df[,-1]) == 1))

  core_otus <- sum(rowSums(membership_df[,-1]) == 4)
  cat("\nCore OTUs shared across all 4 groups:", core_otus, "\n")

  # Rank-abundance (Whittaker) curve across Vegetation types (Veg)
  # 1. For each Veg group, pool OTU counts across all samples in that group,
  # then rank OTUs by total abundance (descending) and compute their
  # relative abundance (%) for the rank-abundance curve
  veg_groups <- levels(sample_data(fungi_physeq)$Veg)

  rank_abund_df <- lapply(veg_groups, function(v) {
      samples_in_group <- sample_names(fungi_physeq)[sample_data(
          fungi_physeq)$Veg == v]
      ps_sub <- prune_samples(samples_in_group, fungi_physeq)
      ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)

      otu_totals <- taxa_sums(ps_sub)
      otu_totals <- sort(otu_totals, decreasing = TRUE)

      data.frame(
        Veg = v,
        Rank = seq_along(otu_totals),
        Abundance = as.numeric(otu_totals),
        RelAbundance = as.numeric(otu_totals) / sum(otu_totals) * 100
      )
  }) %>% dplyr::bind_rows()

  rank_abund_df$Veg <- factor(rank_abund_df$Veg, levels = veg_groups)

  # 2. Plot — OTU rank (x) vs. log10 relative abundance (y), one line per Veg
  veg_colors <- c(SMF = "blue", MF = "#1FBF13", UMF = "orange", CR = "red")

  p_rankabund <- ggplot(rank_abund_df, aes(x = Rank, y = RelAbundance, color =
        Veg)) +
    geom_line(linewidth = 0.9) +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(values = veg_colors) +
    labs(
      x = "OTU rank (most to least abundant)",
      y = "Relative abundance (%, log scale)",
      color = "Vegetation type",
      title = "Rank-abundance (Whittaker) curves across vegetation types",
      subtitle = "Steeper curves indicate greater dominance by fewer taxa;
      flatter curves indicate more even communities"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 12),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9.5, color = "grey30")
    )

  p_rankabund

  ggsave("rank_abundance_curve_by_Veg.pdf", p_rankabund, width = 9, height =
      6.5)
  ggsave("rank_abundance_curve_by_Veg.png", p_rankabund, width = 9, height =
      6.5, dpi = 300)

  # 3. Export the underlying data table (useful for supplementary material
  # or for computing formal evenness metrics later, e.g. Pielou's J)
  write.csv(rank_abund_df, "rank_abundance_data_by_Veg.csv", row.names = FALSE)

  # Combine all graphs on singlplot_df

  # Combine all four figures into a single multi-panel figure with patchwork
  # Combine all four figures into a single multi-panel figure — titles/subtitles
  # removed from each panel (redundant with A/B/C/D tags; captions live in text)

  # 1. Strip titles/subtitles, keep legends compact for multi-panel use
  p_combined_panel <- p_combined +
    labs(title = NULL, subtitle = NULL) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    )

  p_guild_panel <- p_guild +
    labs(title = NULL, subtitle = NULL) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 7.5),
      legend.title = element_text(size = 9)
    )

  p_venn_panel <- p_venn +
    labs(title = NULL)

  p_rankabund_panel <- p_rankabund +
    labs(title = NULL, subtitle = NULL) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    )

  # 2. Assemble as a 2x2 grid with automatic A/B/C/D tags
  final_figure <- (p_combined_panel | p_guild_panel) /
    (p_venn_panel | p_rankabund_panel) +
    plot_annotation(tag_levels = "a", tag_suffix = ")") &
    theme(plot.tag = element_text(face = "bold", size = 15))

  cairo_pdf(
    "combined_multipanel_figure.pdf",
    width = 16,
    height = 12
  )

  print(final_figure)

  dev.off()

  ggsave("combined_multipanel_figure.png", final_figure, width = 16, height =
      12, dpi = 300)

  # Relative abundance stacked bar charts across Veg — Phylum/Order/Family/Genus

  # ----------------------------------------------------------------------------
  # 08. Taxonomic composition
  # ----------------------------------------------------------------------------
  set.seed(123)

  # 1. Convert to relative abundance BEFORE agglomerating at any rank
  fungi_rel <- transform_sample_counts(fungi_physeq, function(x) x / sum(x))

  # 2. Reusable function: build a relative-abundance stacked bar for one rank
  # 0.001 = 0.1%) for a taxon to be shown individually; anything below is
  # pooled into "Other". "Unclassified" is always kept as its own category.
  build_composition_plot <- function(ps_rel, rank, threshold = 0.001,
    veg_levels) {

    other_label <- paste0("Other (<", threshold * 100, "%)")

    glom <- tax_glom(ps_rel, taxrank = rank, NArm = FALSE)
    df <- psmelt(glom)

    df[[rank]] <- as.character(df[[rank]])
    df[[rank]][is.na(df[[rank]]) | df[[rank]] == ""] <- "Unclassified"

    # label (especially "Unclassified", where many different lineages collapse
    # to the same string at this rank). Sum them into ONE value per
    # Sample × rank-label FIRST, so each sample contributes exactly one row —
    # otherwise mean() across samples divides by an inflated row count and
    # silently understates that label's true contribution.
    df_persample <- df %>%
      dplyr::group_by(Sample, Veg, .data[[rank]]) %>%
      dplyr::summarise(Abundance = sum(Abundance), .groups = "drop")

    summary_df <- df_persample %>%
      dplyr::group_by(Veg, .data[[rank]]) %>%
      dplyr::summarise(mean_RA = mean(Abundance), .groups = "drop")

    overall_mean <- summary_df %>%
      dplyr::group_by(.data[[rank]]) %>%
      dplyr::summarise(overall = mean(mean_RA), .groups = "drop")

    keep_taxa <- overall_mean %>%
      dplyr::filter(.data[[rank]] != "Unclassified", overall >= threshold) %>%
      dplyr::pull(.data[[rank]])

    summary_df <- summary_df %>%
      dplyr::mutate(
        TaxonGrouped = case_when(
          .data[[rank]] == "Unclassified" ~ "Unclassified",
          .data[[rank]] %in% keep_taxa ~ .data[[rank]],
          TRUE ~ other_label
        )
      ) %>%
      dplyr::group_by(Veg, TaxonGrouped) %>%
      dplyr::summarise(mean_RA = sum(mean_RA), .groups = "drop")

    abundance_rank <- summary_df %>%
      dplyr::filter(!TaxonGrouped %in% c(other_label, "Unclassified")) %>%
      dplyr::group_by(TaxonGrouped) %>%
      dplyr::summarise(total = sum(mean_RA), .groups = "drop") %>%
      dplyr::arrange(desc(total)) %>%
      dplyr::pull(TaxonGrouped)

    taxon_order <- c(abundance_rank, other_label, "Unclassified")
    taxon_order <- intersect(taxon_order, unique(summary_df$TaxonGrouped))

    stopifnot(all(unique(summary_df$TaxonGrouped) %in% taxon_order))

    summary_df$TaxonGrouped <- factor(summary_df$TaxonGrouped, levels =
        taxon_order)
    summary_df$Veg <- factor(summary_df$Veg, levels = veg_levels)

    totals_check <- summary_df %>% dplyr::group_by(Veg) %>% dplyr::summarise(
      total = sum(mean_RA))
    if (any(abs(totals_check$total - 1) > 0.001)) {
      warning(paste0("Rank '", rank,
          "': group totals do not sum to 1.0 — check for data loss.\n"),
        paste(capture.output(print(totals_check)), collapse = "\n"))
    }

    n_taxa <- length(taxon_order) - sum(c(other_label, "Unclassified") %in%
        taxon_order)
    pal <- as.vector(Polychrome::createPalette(
        max(n_taxa, 3), c("#010D8E", "#F28E2B", "#00A087", "#A65628","#59A14F",
          "#EB0004", "#670A7F"), M = 1000
    ))[seq_len(n_taxa)]
    names(pal) <- taxon_order[seq_len(n_taxa)]
    if (other_label %in% taxon_order)     pal[other_label] <- "grey75"
    if ("Unclassified" %in% taxon_order)  pal["Unclassified"] <- "grey45"

    p <- ggplot(summary_df, aes(x = Veg, y = mean_RA, fill = TaxonGrouped)) +
      geom_col(width = 0.7, color = "white", linewidth = 0.2) +
      scale_fill_manual(values = pal) +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      labs(x = NULL, y = "Mean relative abundance", fill = rank,
        title = paste0(rank, "-level composition (\u2265", threshold * 100,
          "% shown individually)")) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold", size = 10),
        legend.position = "right",
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold", size = 11.5)
      )

    list(plot = p, data = summary_df)
  }
  # Build the five rank-level panels; retain analysis warnings.
  veg_levels <- levels(sample_data(fungi_physeq)$Veg)

  phylum_result <- build_composition_plot(fungi_rel, "Phylum", threshold =
      0.001, veg_levels)
  class_result  <- build_composition_plot(fungi_rel, "Class",  threshold =
      0.001, veg_levels)
  order_result  <- build_composition_plot(fungi_rel, "Order",  threshold =
      0.001, veg_levels)
  family_result <- build_composition_plot(fungi_rel, "Family", threshold =
      0.001, veg_levels)
  genus_result  <- build_composition_plot(fungi_rel, "Genus",  threshold =
      0.001, veg_levels)

  # Combine the five rank-level panels into two rows.
  combined_composition <- (phylum_result$plot | class_result$plot |
      order_result$plot) /
    (family_result$plot | genus_result$plot) +
    plot_annotation(tag_levels = "a", tag_suffix = ")") &
    theme(plot.tag = element_text(face = "bold", size = 14))

  combined_composition

  ggsave(
    "relative_abundance_composition_all_ranks.pdf",
    plot = combined_composition,
    device = cairo_pdf,
    width = 16,
    height = 11,
    units = "in"
  )

  # 5. Export underlying data tables (useful supplementary material)
  write.csv(phylum_result$data, "composition_Phylum_by_Veg.csv", row.names =
      FALSE)
  write.csv(class_result$data,  "composition_Class_by_Veg.csv",  row.names =
      FALSE)
  write.csv(order_result$data,  "composition_Order_by_Veg.csv",  row.names =
      FALSE)
  write.csv(family_result$data, "composition_Family_by_Veg.csv", row.names =
      FALSE)
  write.csv(genus_result$data,  "composition_Genus_by_Veg.csv",  row.names =
      FALSE)

  # 6. Final confirmation — print totals per Veg per rank, should all be ~1.0
  cat("=== Final totals check (should all be ~1.0 / 100%) ===\n")
  list(Phylum = phylum_result$data, Class = class_result$data, Order =
      order_result$data,
    Family = family_result$data, Genus = genus_result$data) %>%
    purrr::imap_dfr(~ .x %>% dplyr::group_by(Veg) %>% dplyr::summarise(total =
          sum(mean_RA)) %>% dplyr::mutate(Rank = .y)) %>%

  # ----------------------------------------------------------------------------
  # 09. Alpha diversity and Hill profiles
  # ----------------------------------------------------------------------------
    print(n = Inf)
  # Hill numbers and community ordination

  set.seed(123)
  # )
  veg_colors <- c(SMF = "orange", MF = "blue", UMF = "green3", CR = "red" )
  veg_levels <- levels(sample_data(fungi_physeq)$Veg)

  # PART 1 — Hill numbers (q = 0, 1, 2)

  # 1. Extract OTU table as samples x taxa proportions (Hill numbers use
  # proportions, NOT CLR-transformed data — that's for the beta-diversity
  # section below)
  otu_mat <- as(otu_table(fungi_physeq), "matrix")
  # samples x taxa
  if (taxa_are_rows(fungi_physeq)) otu_mat <- t(otu_mat)

  prop_mat <- otu_mat / rowSums(otu_mat)

  # 2. Compute Hill numbers per sample
  hill_q0 <- rowSums(prop_mat > 0)

  hill_q1 <- apply(prop_mat, 1, function(p) {
      p <- p[p > 0]
      exp(-sum(p * log(p)))
  })

  hill_q2 <- apply(prop_mat, 1, function(p) {
      1 / sum(p^2)
  })

  hill_df <- data.frame(
    Sample = rownames(prop_mat),
    q0 = hill_q0,
    q1 = hill_q1,
    q2 = hill_q2
  ) %>%
    dplyr::left_join(
      as(sample_data(fungi_physeq), "data.frame") %>%
        tibble::rownames_to_column("Sample") %>%
        dplyr::select(Sample, Veg),
      by = "Sample"
    ) %>%
    dplyr::mutate(Veg = factor(Veg, levels = veg_levels))

  write.csv(hill_df, "hill_numbers_by_sample.csv", row.names = FALSE)

  # 3. Reshape long, run Kruskal-Wallis + Dunn's post-hoc (FDR-corrected)
  # per Hill order, then build a faceted boxplot with significance brackets
  hill_long <- hill_df %>%
    pivot_longer(cols = c(q0, q1, q2), names_to = "HillOrder", values_to =
        "Value") %>%
    dplyr::mutate(HillOrder = factor(HillOrder, levels = c("q0", "q1", "q2"),
        labels = c("q=0 (Richness)", "q=1 (Shannon)", "q=2 (Simpson)")))

  kw_results <- hill_long %>%
    dplyr::group_by(HillOrder) %>%
    kruskal_test(Value ~ Veg)

  cat("=== Kruskal-Wallis results across Veg, per Hill order ===\n")
  print(kw_results)

  dunn_results_clean <- hill_long %>%
    dplyr::group_by(HillOrder) %>%
    dunn_test(Value ~ Veg, p.adjust.method = "BH")

  write.csv(dunn_results_clean, "hill_dunn_posthoc_results.csv", row.names =
      FALSE)

  # Separate copy WITH plot-position columns, used only for stat_pvalue_manual()

  # Improved faceted Hill number plot — tighter bracket placement so boxplots
  # aren't squished by y-axis inflation from significance annotations

  # 1. Recompute bracket positions with a SMALLER step increase (was 0.08),
  # and base the starting bracket height tighter to the actual data max
  # per facet, rather than the default generous spacing
  dunn_results <- dunn_results_clean %>%
    add_xy_position(x = "Veg", step.increase = 0.04, scales = "free_y")

  # 2. Plot — tighter y-axis expansion so the top of the panel doesn't leave
  # a large empty gap above the highest bracket
  p_hill <- ggplot(hill_long, aes(x = Veg, y = Value, fill = Veg)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.55) +
    geom_jitter(width = 0.12, size = 3.5, alpha = 0.6, stroke = 0.3, shape =
        21) +
    facet_wrap(~ HillOrder, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = veg_colors) +
    stat_pvalue_manual(dunn_results, hide.ns = TRUE, tip.length = 0.10,
      size = 7, bracket.size = 0.4) +
  # tighter top margin
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
    labs(x = "Vegetation type", y = "Hill number value", fill = NULL) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey95", color = "grey30",
        linewidth = 0.6),
      axis.text.x = element_text(face = "bold", size = 10),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14)
    )

  p_hill

  ggsave(
    "hill_numbers_by_Veg_v2.pdf",
    plot = p_hill,
    device = cairo_pdf,
    width = 16,
    height = 11,
    units = "in"
  )
  # plot — half-violin (distribution) + boxplot (summary) +
  # jittered points (raw data), replacing the plain boxplot

  p_hill_sina <- ggplot(hill_long, aes(x = Veg, y = Value, fill = Veg)) +
    geom_violin(alpha = 0.3, color = NA) +
    ggforce::geom_sina(aes(color = Veg), size = 2.5, alpha = 0.8, maxwidth =
        0.6) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, color =
        "black") +
    facet_wrap(~ HillOrder, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = veg_colors) +
    scale_color_manual(values = veg_colors) +
    stat_pvalue_manual(dunn_results, hide.ns = TRUE, tip.length = 0.01, size =
        6, bracket.size = 0.4) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = "Vegetation type", y = "Hill number value") +
    theme_bw(base_size = 13) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(),
      strip.text = element_text(face = "bold"), plot.title = element_text(face =
          "bold", size = 14))

  p_hill_sina

  # Hill number diversity PROFILE — continuous q, one curve per sample
  # 1. Compute Hill numbers across a continuous range of q (0 to 2, fine steps)
  q_seq <- seq(0, 2, by = 0.1)

  compute_hill_q <- function(p, q) {
    p <- p[p > 0]
    if (q == 1) {
      exp(-sum(p * log(p)))
    } else {
      (sum(p^q))^(1 / (1 - q))
    }
  }

  hill_profile <- expand.grid(Sample = rownames(prop_mat), q = q_seq) %>%
    rowwise() %>%
    dplyr::mutate(Hill = compute_hill_q(prop_mat[Sample, ], q)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
      as(sample_data(fungi_physeq), "data.frame") %>%
        tibble::rownames_to_column("Sample") %>% dplyr::select(Sample, Veg),
      by = "Sample"
    ) %>%
    dplyr::mutate(Veg = factor(Veg, levels = veg_levels))

  # 2. Group-level mean profile + confidence ribbon (mean ± SE per q, per group)
  hill_profile_summary <- hill_profile %>%
    dplyr::group_by(Veg, q) %>%
    dplyr::summarise(mean_Hill = mean(Hill), se_Hill = sd(Hill) / sqrt(n()),
      .groups = "drop")

  # 3. Plot — thin individual sample curves (background) + bold mean curve
  # with ribbon (foreground), faceted or all-in-one, per group color
  p_hill_profile <- ggplot() +
    geom_line(data = hill_profile, aes(x = q, y = Hill, group = Sample, color =
          Veg),
      alpha = 0.25, linewidth = 0.4) +
    geom_ribbon(data = hill_profile_summary,
      aes(x = q, ymin = mean_Hill - se_Hill, ymax = mean_Hill + se_Hill, fill =
          Veg),
      alpha = 0.25, color = NA) +
    geom_line(data = hill_profile_summary, aes(x = q, y = mean_Hill, color =
          Veg),
      linewidth = 1.3) +
    scale_color_manual(values = veg_colors) +
    scale_fill_manual(values = veg_colors) +
    labs(
      x = "Diversity order (q)", y = "Hill number",
      color = "Vegetation type", fill = "Vegetation type"
    ) +
    annotate("text", x = 0, y = Inf, label = "q=0: Richness", hjust = -0.05,
      vjust = 1.5, size = 3, fontface = "italic", color = "grey40") +
    annotate("text", x = 2, y = Inf, label = "q=2: Dominance-weighted", hjust =
        1.05, vjust = 1.5, size = 3, fontface = "italic", color = "grey40") +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      panel.grid.minor = element_blank()
    )

  p_hill_profile <- p_hill_profile +
    scale_y_log10() +
    labs(y = "Hill number (log scale)")

  p_hill_profile

  ggsave("hill_diversity_profile_log.pdf", p_hill_profile, device = cairo_pdf,
    width = 10, height = 7, units = "in")
  ggsave("hill_diversity_profile_log.png", p_hill_profile, width = 10, height =
      7, dpi = 300, bg = "white")

  # ----------------------------------------------------------------------------
  # 10. Ordination and environmental models
  # ----------------------------------------------------------------------------
  # Rebuild the community matrix at OTU level. This deliberately replaces
  # the family-level clr_mat used only in the earlier correlation heatmap.
  otu_ordination <- as(otu_table(fungi_physeq), "matrix")
  if (taxa_are_rows(fungi_physeq)) {
    otu_ordination <- t(otu_ordination)
  }
  sample_order <- sample_names(fungi_physeq)
  otu_ordination <- otu_ordination[sample_order, , drop = FALSE]
  otu_ordination <- otu_ordination[
    , colSums(otu_ordination) > 0, drop = FALSE
  ]
  if (any(!is.finite(otu_ordination)) || any(otu_ordination < 0)) {
    stop("The OTU matrix must contain finite, non-negative counts.")
  }
  if (nrow(otu_ordination) < 6L || ncol(otu_ordination) < 3L) {
    stop("Ordination requires at least six samples and three detected OTUs.")
  }
  if (!is.finite(clr_pseudocount) || clr_pseudocount <= 0) {
    stop("clr_pseudocount must be a positive finite number.")
  }

  log_abundance <- log(otu_ordination + clr_pseudocount)
  clr_mat <- sweep(log_abundance, 1, rowMeans(log_abundance), "-")
  aitchison_dist <- stats::dist(clr_mat, method = "euclidean")
  hellinger_mat <- vegan::decostand(otu_ordination, method = "hellinger")
  bray_dist <- vegan::vegdist(hellinger_mat, method = "bray")

  # Align environmental data explicitly to the community-matrix sample order.
  metadata_ordination <- as(sample_data(fungi_physeq), "data.frame")
  metadata_ordination <- metadata_ordination[sample_order, , drop = FALSE]
  env_all <- metadata_ordination %>%
    dplyr::select(where(is.numeric)) %>%
    dplyr::select(-any_of(c("Latitude", "Longitude", "MO")))
  finite_variance <- vapply(env_all, function(x) {
      x <- x[is.finite(x)]
      length(x) > 1L && stats::var(x) > 1e-6
    }, logical(1))
  env_all <- env_all[, finite_variance, drop = FALSE]
  if (!ncol(env_all)) stop("No variable environmental predictors remain.")
  complete_env <- apply(as.matrix(env_all), 1, function(x) {
      all(is.finite(x))
  })
  env_matrix <- env_all[complete_env, , drop = FALSE]
  if (nrow(env_matrix) <= ncol(env_matrix) + 2L) {
    stop("Too few complete samples for the full environmental RDA model.")
  }
  write.csv(
    data.frame(Sample = rownames(env_all), Included = complete_env),
    "ordination_environment_sample_inclusion.csv", row.names = FALSE
  )

  # PCoA retains every sample with valid fungal read counts.
  pcoa_result <- stats::cmdscale(aitchison_dist, k = 2, eig = TRUE)
  if (ncol(pcoa_result$points) < 2L) {
    stop("The community data support fewer than two PCoA axes.")
  }
  pcoa_points <- pcoa_result$points[, 1:2, drop = FALSE]
  colnames(pcoa_points) <- c("PCo1", "PCo2")
  positive_eigenvalues <- pcoa_result$eig[pcoa_result$eig > 0]
  pcoa_var_explained <- round(
    100 * pcoa_result$eig[1:2] / sum(positive_eigenvalues), 2
  )
  pcoa_df <- data.frame(
    Sample = rownames(pcoa_points), pcoa_points,
    Veg = metadata_ordination[rownames(pcoa_points), "Veg"]
  )

  set.seed(random_seed)
  nmds_result <- vegan::metaMDS(
    bray_dist, k = 2, trymax = nmds_trymax,
    autotransform = FALSE, trace = FALSE
  )
  nmds_points <- vegan::scores(nmds_result, display = "sites", choices = 1:2)
  colnames(nmds_points) <- c("NMDS1", "NMDS2")
  nmds_df <- data.frame(
    Sample = rownames(nmds_points), nmds_points,
    Veg = metadata_ordination[rownames(nmds_points), "Veg"]
  )
  if (isTRUE(nmds_result$converged == 0L)) {
    warning("NMDS did not repeat its best solution; inspect convergence.")
  }

  # Fit environmental vectors to the two displayed PCoA axes.
  set.seed(random_seed)
  envfit_pcoa <- vegan::envfit(
    pcoa_points[rownames(env_matrix), , drop = FALSE],
    env_matrix, permutations = ordination_permutations
  )
  envfit_vectors <- vegan::scores(envfit_pcoa, display = "vectors")
  envfit_pcoa_df <- data.frame(
    Variable = rownames(envfit_vectors), envfit_vectors,
    r2 = unname(envfit_pcoa$vectors$r[rownames(envfit_vectors)]),
    p = unname(envfit_pcoa$vectors$pvals[rownames(envfit_vectors)])
  )
  write.csv(pcoa_df, "PCoA_sample_scores.csv", row.names = FALSE)
  write.csv(nmds_df, "NMDS_sample_scores.csv", row.names = FALSE)
  write.csv(envfit_pcoa_df, "PCoA_environment_envfit.csv", row.names = FALSE)
  write.csv(
    data.frame(Axis = c("PCo1", "PCo2"), Percent = pcoa_var_explained),
    "PCoA_variance_explained.csv", row.names = FALSE
  )
  write.csv(
    data.frame(Stress = nmds_result$stress),
    "NMDS_stress.csv", row.names = FALSE
  )
  # This is an output cache, not an additional input requirement.
  saveRDS(
    list(
      clr_mat = clr_mat, env_matrix = env_matrix,
      aitchison_dist = aitchison_dist, bray_dist = bray_dist,
      pcoa_result = pcoa_result, pcoa_df = pcoa_df,
      pcoa_var_explained = pcoa_var_explained,
      nmds_result = nmds_result, nmds_df = nmds_df,
      envfit_pcoa_df = envfit_pcoa_df
    ),
    "computed_ordination_objects.rds"
  )

  # VIF screening for the primary environmental RDA model.
  env_complete <- env_matrix[complete.cases(env_matrix), , drop = FALSE]
  clr_complete <- clr_mat[rownames(env_complete), , drop = FALSE]
  cat("\nSamples retained for RDA:", nrow(env_complete), "of", nrow(
      env_matrix), "\n")

  # DIAGNOSTIC ONLY — never reported
  rda_full <- rda(clr_complete ~ ., data = env_complete)
  cat("\n=== Initial VIF (full 17-variable model — diagnostic step) ===\n")
  print(sort(vif.cca(rda_full), decreasing = TRUE))

  remove_high_vif <- function(env_data, response_data, threshold = 10) {
    removed_log <- tibble(Variable = character(), VIF = numeric())
    repeat {
      mod <- rda(response_data ~ ., data = env_data)
      v <- vif.cca(mod)
      v[!is.finite(v)] <- Inf
      if (length(v) <= 1 || all(v < threshold)) break
      worst <- names(which.max(v))
      cat("Removing (VIF =", round(max(v), 2), "):", worst, "\n")
      removed_log <- dplyr::bind_rows(removed_log, tibble(Variable = worst,
          VIF = max(v)))
      env_data <- env_data[, setdiff(colnames(env_data), worst), drop = FALSE]
    }
    list(env_data = env_data, model = rda(response_data ~ ., data = env_data),
      removed = removed_log)
  }

  vif_reduced <- remove_high_vif(env_complete, clr_complete, threshold = 10)
  env_reduced <- vif_reduced$env_data

  cat(paste0("\n=== Variables retained after VIF ",
      "reduction (13 variables — PRIMARY model) ", "===\n"))
  print(colnames(env_reduced))
  write.csv(vif_reduced$removed, "SuppTable_VIF_removed_variables.csv",
    row.names = FALSE)

  # Fit the primary RDA using all retained environmental constraints.
  rda_model_reduced <- rda(clr_complete ~ ., data = env_reduced)

  cat("\n=== Primary RDA model with", ncol(env_reduced), "predictors ===\n")
  print(rda_model_reduced)

  rda_r2_final <- RsquareAdj(rda_model_reduced)
  rda_anova_final <- anova.cca(rda_model_reduced, permutations = 999)
  rda_axis_anova_final <- anova.cca(rda_model_reduced, by = "axis",
    permutations = 999)
  rda_term_anova_final <- anova.cca(rda_model_reduced, by = "margin",
    permutations = 999)

  cat("\n=== Overall RDA significance ===\n"); print(rda_anova_final)
  cat("\n=== Per-axis RDA significance ===\n"); print(rda_axis_anova_final)
  cat(paste0("\n=== Per-variable (marginal) ",
      "significance — Supplementary Table ===\n")); print(rda_term_anova_final)
  cat("\nAdjusted R\u00b2:", round(rda_r2_final$adj.r.squared, 4), "\n")

  table_rda_marginal <- as.data.frame(rda_term_anova_final) %>%
    rownames_to_column("Variable") %>%
    dplyr::filter(Variable != "Residual") %>%
    dplyr::transmute(Variable, Variance = round(Variance, 1), F = round(F, 3),
      `p-value` = `Pr(>F)`) %>%
    dplyr::arrange(desc(Variance))

  sig_vars_final <- table_rda_marginal$Variable[table_rda_marginal$`p-value` <
    0.05]
  cat("\nSignificant variables:", paste(sig_vars_final, collapse = ", "),
    "\n")

  write.csv(table_rda_marginal,
    "SuppTable_RDA_13variable_marginal_significance.csv", row.names = FALSE)

  # 4. Adjusted R² as its own small summary file (since it's a single value,
  # not a data frame like the other three)
  r2_summary <- data.frame(
    Metric = c("R2", "Adjusted R2"),
    Value = c(round(rda_r2_final$r.squared, 4), round(
        rda_r2_final$adj.r.squared, 4))
  )
  write.csv(r2_summary, "RDA_R2_summary.csv", row.names = FALSE)
  write.csv(as.data.frame(rda_anova_final),
    "RDA_overall_significance.csv")
  write.csv(as.data.frame(rda_axis_anova_final),
    "RDA_per_axis_significance.csv")
  write.csv(as.data.frame(rda_term_anova_final),
    "RDA_per_variable_marginal_significance.csv")
  sink("RDA_full_significance_summary.csv")
  cat("=== Overall RDA significance ===\n")
  write.csv(as.data.frame(rda_anova_final), row.names = TRUE)
  cat("\n=== Per-axis RDA significance ===\n")
  write.csv(as.data.frame(rda_axis_anova_final), row.names = TRUE)
  cat("\n=== Per-variable (marginal) significance ===\n")
  write.csv(as.data.frame(rda_term_anova_final), row.names = TRUE)
  cat("\n=== R-squared summary ===\n")
  write.csv(r2_summary, row.names = FALSE)
  sink()

  cat(paste0("Saved: RDA_overall_significance.csv, ",
      "RDA_per_axis_significance.csv,\n"))
  cat(paste0("       RDA_per_variable_marginal_significa",
      "nce.csv, RDA_R2_summary.csv,\n"))
  cat("       RDA_full_significance_summary.csv (combined)\n")

  # Extract site scores and the retained environmental biplot vectors.
  veg_lookup <- as(sample_data(fungi_physeq), "data.frame") %>%
    rownames_to_column("Sample_orig")
  veg_for_rda <- veg_lookup$Veg[match(rownames(env_reduced),
      veg_lookup$Sample_orig)]

  rda_scores_sites <- as.data.frame(scores(rda_model_reduced, display =
        "sites")) %>%
    rownames_to_column("Sample") %>%
    dplyr::mutate(Veg = factor(veg_for_rda, levels = veg_levels))

  rda_scores_env_final <- as.data.frame(scores(rda_model_reduced, display =
        "bp")) %>%
    rownames_to_column("Variable")

  # Sensitivity comparison with a forward-selected model.
  # Axis agreement alone does not establish absence of overfitting.
  rda_null <- rda(clr_complete ~ 1, data = env_reduced)
  rda_selected <- ordiR2step(rda_null, scope = formula(rda_model_reduced),
    direction = "forward", permutations = 999, trace = TRUE)

  cat("\n=== Robustness check: forward-selected model ===\n")
  print(rda_selected$call)

  robustness_cor <- NA_real_
  robustness_maxdiff <- NA_real_
  if (!is.null(rda_selected$CCA) && rda_selected$CCA$rank > 0L) {
    scores_full <- vegan::scores(rda_model_reduced, display = "sites")
    scores_selected <- vegan::scores(rda_selected, display = "sites")
    robustness_cor <- cor(scores_full[, 1], scores_selected[, 1])
    robustness_maxdiff <- max(abs(scores_full[, 1] - scores_selected[, 1]))
  } else {
    warning("Forward selection retained no constraints; comparison is NA.")
  }

  cat("\nCorrelation between 13-variable and 3-variable RDA1 scores:", round(
      robustness_cor, 4), "\n")
  cat("Max absolute difference:", round(robustness_maxdiff, 4), "\n")

  robustness_table <- tibble(
    Comparison =
      paste0(ncol(env_reduced), "-variable VIF model vs. ",
        if (is.null(rda_selected$CCA)) 0L else rda_selected$CCA$qrank,
        "-variable selected model"),
    `Pearson's r (RDA1)` = round(robustness_cor, 4),
    `Max absolute difference` = round(robustness_maxdiff, 4)
  )
  write.csv(robustness_table, "SuppTable_RDA_model_robustness_comparison.csv",
    row.names = FALSE)
  # required colors for the plot
  veg_colors <- c(SMF = "orange", MF = "blue", UMF = "green3", CR = "red")
  veg_levels <- levels(sample_data(fungi_physeq)$Veg)

  # SECTION G — MAIN-TEXT FIGURE: PCoA (a) | NMDS (b) | RDA biplot (c)
  # Ellipses filled + points filled, both colored by veg_colors
  p_pcoa <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2)) +
    stat_ellipse(aes(color = Veg, fill = Veg), geom = "polygon", alpha = 0.20,
      level = 0.95, linewidth = 0.7) +
    geom_point(aes(fill = Veg), shape = 21, color = "black", size = 3,
      stroke = 0.4) +
    scale_color_manual(values = veg_colors) +
    scale_fill_manual(values = veg_colors) +
    labs(title = "PCoA (Aitchison distance)",
      x = paste0("PCo1 (", pcoa_var_explained[1], "%)"), y = paste0("PCo2 (",
        pcoa_var_explained[2], "%)")) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12.5),
      legend.position = "none",
      panel.grid.minor = element_blank())

  p_nmds <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
    stat_ellipse(aes(color = Veg, fill = Veg), geom = "polygon", alpha = 0.20,
      level = 0.95, linewidth = 0.7) +
    geom_point(aes(fill = Veg), shape = 21, color = "black", size = 3,
      stroke = 0.4) +
    scale_color_manual(values = veg_colors) +
    scale_fill_manual(values = veg_colors) +
    labs(title = paste0("NMDS, Bray-Curtis (stress = ", round(
    nmds_result$stress, 3), ")"),
      x = "NMDS1", y = "NMDS2") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12.5),
      legend.position = "none",
      panel.grid.minor = element_blank())

  arrow_scale <- 9
  p_rda <- ggplot(rda_scores_sites, aes(x = RDA1, y = RDA2)) +
    stat_ellipse(aes(color = Veg, fill = Veg), geom = "polygon", alpha = 0.20,
      level = 0.95, linewidth = 0.7) +
    geom_point(aes(fill = Veg), shape = 21, color = "black", size = 3,
      stroke = 0.4)+
    geom_segment(data = rda_scores_env_final,
      aes(x = 0, y = 0, xend = RDA1 * arrow_scale, yend = RDA2 * arrow_scale),
      inherit.aes = FALSE, arrow = arrow(length = unit(0.15, "cm")), color =
        "grey20", linewidth = 0.5) +
    ggrepel::geom_text_repel(data = rda_scores_env_final,
      aes(x = RDA1 * arrow_scale, y = RDA2 * arrow_scale, label = Variable),
      inherit.aes = FALSE, size = 3.4, color = "grey20", max.overlaps = 20) +
    scale_color_manual(values = veg_colors) +
    scale_fill_manual(values = veg_colors) +
    labs(title = paste0("RDA (Adj. R\u00b2 = ", round(
          rda_r2_final$adj.r.squared, 3), ")"),
      subtitle = paste0("Constrained by ", ncol(env_reduced),
        " predictors (VIF < 10)"),
      x = "RDA1", y = "RDA2") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12.5),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      panel.grid.minor = element_blank())

  combined_ordination <- p_pcoa + p_nmds + p_rda +
    plot_layout(guides = "collect", nrow = 1) +
    plot_annotation(title = "Community ordination and environmental drivers",
      tag_levels = "a", tag_suffix = ")") &
    theme(legend.position = "right", plot.tag = element_text(face = "bold",
        size = 13),
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5))

  combined_ordination
  ggsave("Figure_ordination_PCoA_NMDS_RDA.pdf", combined_ordination, device =
      cairo_pdf, width = 18, height = 6.5)
  ggsave("Figure_ordination_PCoA_NMDS_RDA.png", combined_ordination, width =
      18, height = 6.5, dpi = 300, bg = "white")
  # Combine the p-hill plot with combined_ordination_final using patchwork
  diversity_and_ordination <- (p_hill_sina + p_hill_profile) /
    combined_ordination +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(
      tag_levels = "a", tag_suffix = ")"
    ) &
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.tag = element_text(face = "bold", size = 13)
    )
  diversity_and_ordination
  # Save diversity_and_ordination as cairo pdf

  ggsave("Figure_diversity_and_ordination.pdf", diversity_and_ordination,
    device = cairo_pdf, width = 18, height = 13, units = "in")

  # ----------------------------------------------------------------------------
  # 11. Genus association networks and visualisation
  # ----------------------------------------------------------------------------
  # Genus-level, full taxonomic rank + guild attachment, circle-packing plots

  set.seed(123)

  veg_levels <- levels(sample_data(fungi_physeq)$Veg)

  # "Mixed" = "#8C6D31", "Generalist saprotroph" = "#EF8A17", "Endophyte" =
  # "#519E3E",
  # "Plant pathogen" = "#C53A32", "Wood saprotroph" = "#8E5B9E", "Dung
  # saprotroph" = "#B07AA1",
  # "Plant saprotroph" = "#D4B02C", "Mycorrhizal" = "#1B9E77", "Fungal parasite"
  # = "#E377C2",
  # "Animal pathogen" = "#D6616B", "Lichen parasite" = "#17BECF", "Lichenized" =
  # "#B5CF6B",
  # "Nematophagous" = "#7B4173", "Generalist pathogen" = "#843C39", "Unassigned"
  # = "grey70"
  # )

  # SECTION A — Build genus-level table + guild lookup (shared across all
  # groups)
  fungi_genus <- tax_glom(fungi_physeq, taxrank = "Genus", NArm = FALSE)
  genus_labels <- as.data.frame(tax_table(fungi_genus))$Genus
  taxa_names(fungi_genus) <- make.unique(as.character(genus_labels))

  # tax_glom does not reliably carry custom columns (like Guild_class) through
  # agglomeration, so build a genus -> guild lookup directly from the
  # OTU-level fungi_physeq tax_table instead of trusting fungi_genus's copy
  otu_tax_df <- as.data.frame(tax_table(fungi_physeq))
  guild_lookup <- otu_tax_df %>%
    dplyr::filter(!is.na(Genus), Genus != "Unclassified") %>%
    dplyr::distinct(Genus, .keep_all = TRUE) %>%
    dplyr::select(Genus, Guild_class)

  # SECTION B — Single build function, one network per Veg group
  build_veg_network <- function(veg_group,
    prevalence_thresh = 0.50,
    min_total_count = 20,
    lambda.min.ratio = 1e-3,
    nlambda = 30,
    rep.num = 50) {

    cat("\n========== Building network for:", veg_group, "==========\n")

    samples_in_group <- sample_names(fungi_genus)[sample_data(
        fungi_genus)$Veg == veg_group]
    ps_sub <- prune_samples(samples_in_group, fungi_genus)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)

    otu_mat <- as(otu_table(ps_sub), "matrix")
    if (!taxa_are_rows(ps_sub)) otu_mat <- t(otu_mat)

    prevalence <- rowSums(otu_mat > 0) / ncol(otu_mat)
    total_counts <- rowSums(otu_mat)

    keep_taxa <- names(prevalence)[prevalence >= prevalence_thresh &
        total_counts >= min_total_count]
    ps_filt <- prune_taxa(keep_taxa, ps_sub)

    cat("Samples:", nsamples(ps_filt), "| Taxa retained:", ntaxa(ps_filt), "\n")

    if (ntaxa(ps_filt) < 5) {
      warning(veg_group,
        ": fewer than 5 taxa survived filtering — network not meaningful.")
      return(NULL)
    }

    se_mb <- SpiecEasi::spiec.easi(
      ps_filt, method = "mb", lambda.min.ratio = lambda.min.ratio, nlambda =
        nlambda,
      pulsar.params = list(rep.num = rep.num, ncores = 1, thresh = 0.05)
    )

    opt_index <- SpiecEasi::getOptInd(se_mb)
    stability <- SpiecEasi::getStability(se_mb)
    cat("StARS stability:", round(stability, 4), "| Selected lambda index:",
      opt_index,
      "of", length(se_mb$est$lambda), "\n")

    if (opt_index %in% c(1, length(se_mb$est$lambda))) {
      warning(veg_group, ": opt_index at grid boundary (", opt_index,
        ") — network may not be fully converged. Consider tightening the ",
        "taxon filter or widening lambda.min.ratio.")
    }

    ig <- SpiecEasi::adj2igraph(SpiecEasi::getRefit(se_mb),
      vertex.attr = list(name = taxa_names(ps_filt)))

    se_beta <- SpiecEasi::symBeta(SpiecEasi::getOptBeta(se_mb), mode = "maxabs")
    edge_weights <- as.matrix(se_beta)
    rownames(edge_weights) <- taxa_names(ps_filt)
    colnames(edge_weights) <- taxa_names(ps_filt)

    E(ig)$weight <- edge_weights[cbind(
        match(igraph::head_of(ig, E(ig))$name, rownames(edge_weights)),
        match(igraph::tail_of(ig, E(ig))$name, colnames(edge_weights))
    )]
    E(ig)$sign <- ifelse(E(ig)$weight > 0, "positive", "negative")

    isolated <- which(igraph::degree(ig) == 0)
    ig <- igraph::delete_vertices(ig, isolated)

    # Guild via direct genus-name join against the verified OTU-level lookup
    base_genus <- sub("\\.[0-9]+$", "", V(ig)$name)
    matched_guild <- guild_lookup$Guild_class[match(base_genus,
        guild_lookup$Genus)]
    V(ig)$Guild_class <- ifelse(is.na(matched_guild) | matched_guild %in% c(
        "Unknown", ""),
      "Unassigned", matched_guild)

    # Attach all taxonomic ranks (Phylum through Species) as node attributes —
    # needed for plot_circlepack()'s flexible group_attr/color_attr options
    tax_df <- as.data.frame(tax_table(ps_filt))
    for (rnk in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
      if (rnk %in% colnames(tax_df)) {
        vertex_attr(ig, rnk) <- tax_df[V(ig)$name, rnk]
      }
    }

    V(ig)$degree <- degree(ig)

    cat("Final network:", vcount(ig), "nodes,", ecount(ig),
      "edges | Guild match rate:",
      sum(V(ig)$Guild_class != "Unassigned"), "of", vcount(ig), "\n")

    list(graph = ig, model = se_mb, ps_filt = ps_filt, stability = stability,
      opt_index = opt_index)
  }

  # Run the build across all four groups
  veg_networks <- purrr::map(veg_levels, build_veg_network)
  names(veg_networks) <- veg_levels
  veg_networks <- veg_networks[!sapply(veg_networks, is.null)]
  if (!all(c("SMF", "MF", "UMF", "CR") %in% names(veg_networks))) {
    stop("A habitat network failed filtering; review before plotting.")
  }
  if (any(vapply(veg_networks, function(x) vcount(x$graph) == 0,
        logical(1)))) {
    stop("A habitat network has no connected vertices after inference.")
  }

  cat("\nGroups with valid networks:", paste(names(veg_networks), collapse =
        ", "), "\n")

  cat("\nColumns present in SMF's node attributes:\n")
  print(colnames(igraph::as_data_frame(veg_networks$SMF$graph, what =
          "vertices")))

  # Circular hierarchical edge bundling — manual guild/phylum palettes,
  # per-group top-N labeling, uniform node size
  # High-contrast manual palettes — mix of dark and bright hues, chosen so
  # no two adjacent legend entries are easily confused
  guild_seed_colors <- c(
    # near-black
    "Mixed"                  = "#222222",
    # bright orange
    "Generalist saprotroph"  = "#F38400",
    # bright green
    "Endophyte"              = "green",
    # bright red
    "Plant pathogen"         = "#BE0032",
    # medium purple
    "Wood saprotroph"        = "purple",
    # dark brown
    "Dung saprotroph"        = "#7F080C",
    # bright yellow
    "Plant saprotroph"       = "#F3C300",
    # bright blue
    "Mycorrhizal"            = "skyblue",
    # bright pink
    "Fungal parasite"        = "pink",
    # dark reddish-brown
    "Animal pathogen"        = "#B83D1F",
    # bright yellow-green
    "Lichen parasite"        = "#DCD300",
    # dark olive green
    "Lichenized"             = "darkgreen",
    # reddish magenta
    "Nematophagous"          = "hotpink",
    # navy (dark)
    "Generalist pathogen"    = "#000080",
    # neutral grey
    "Unassigned"             = "grey"
  )

  phylum_seed_colors <- c(
    # bright red
    "Ascomycota"          = "#BE0032",
    # bright green
    "Basidiomycota"       = "darkgreen",
    # bright blue
    "Mortierellomycota"   = "#377EB8",
    # bright orange
    "Chytridiomycota"     = "orange",
    # bright yellow (colorblind-safe)
    "Glomeromycota"       = "yellow",
    # dark brown
    "Mucoromycota"        = "#7F080C",
    # bright pink
    "Rozellomycota"       = "pink",
    # black
    "Zoopagomycota"       = "#000000",
    # bright teal
    "Basidiobolomycota"   = "green",
    # dark magenta
    "Kickxellomycota"     = "purple",
    # indigo (dark)
    "Monoblepharomycota"  = "#4B0082",
    # dark navy
    "Aphelidiomycota"     = "#08306B",
    # neutral grey
    "Unclassified"        = "#B0B0B0"
  )

  # Diagnostic: confirm palettes cover every level actually present
  check_palette_coverage <- function(veg_networks, group_attr = "Guild_class",
    color_attr = "Phylum") {
    all_guilds <- unique(unlist(lapply(veg_networks, function(n) {
            igraph::vertex_attr(n$graph, group_attr)
    })))
    all_phyla <- unique(unlist(lapply(veg_networks, function(n) {
            igraph::vertex_attr(n$graph, color_attr)
    })))

    missing_guilds <- setdiff(all_guilds, names(guild_seed_colors))
    missing_phyla  <- setdiff(all_phyla, names(phylum_seed_colors))

    if (length(missing_guilds) > 0) cat("WARNING: Guild levels NOT in palette:",
      paste(missing_guilds, collapse = ", "), "\n")
    if (length(missing_phyla)  > 0) cat(
      "WARNING: Phylum levels NOT in palette:", paste(missing_phyla, collapse =
          ", "), "\n")
    if (length(missing_guilds) == 0 && length(missing_phyla) == 0) cat(
      "Palette coverage OK — all levels present are covered.\n")
  }

  check_palette_coverage(veg_networks)

  # separate contiguous segments (a guild can legitimately appear in
  # multiple, non-adjacent clusters around the ring)
  make_palette <- function(values, seed_colors = NULL) {
    levels_present <- sort(unique(na.omit(values)))
    n <- length(levels_present)
    if (!is.null(seed_colors) && all(levels_present %in% names(seed_colors))) {
      return(seed_colors[levels_present])
    }
    pal <- if (n <= 12) scales::hue_pal()(n) else colorRampPalette(
      scales::hue_pal()(12))(n)
    setNames(pal, levels_present)
  }

  plot_circular_bundle <- function(net_obj, veg_group,
    group_attr = "Guild_class",
    color_attr = "Phylum",
    size_attr  = "degree",
    min_group_count = 1,
    min_color_count = 1,
    label_rank = NULL,
    label_mode = "per_group_top_n",
    label_top_n = 10,
    label_size = 3.4,
    exclude_unclassified_labels = TRUE,
    group_colors = NULL,
    node_colors  = NULL,
    show_node_legend = TRUE,
    edge_width = 0.2,
    edge_alpha = 0.3,
    tension = 0.85,
    node_size = 2.5,
    band_inner_r = 1.02,
    band_outer_r = 1.10,
    band_alpha = 0.85) {

    ig <- net_obj$graph
    node_df <- igraph::as_data_frame(ig, what = "vertices")

    stopifnot(all(c(group_attr, color_attr, size_attr) %in% colnames(node_df)))

    node_df$.group <- as.character(node_df[[group_attr]])
    node_df$.group[is.na(node_df$.group) | node_df$.group == ""] <- "Unassigned"

    node_df$.color <- as.character(node_df[[color_attr]])
    node_df$.color[is.na(node_df$.color) | node_df$.color == ""] <-
      "Unclassified"

    # Merge rare categories into "Other" so legends stay meaningful
    group_counts <- table(node_df$.group)
    rare_groups <- names(group_counts)[group_counts < min_group_count]
    if (length(rare_groups) > 0) {
      node_df$.group[node_df$.group %in% rare_groups] <- "Other"
      cat("Merged into 'Other' (guild, n<", min_group_count, "):", paste(
          rare_groups, collapse = ", "), "\n")
    }

    color_counts <- table(node_df$.color)
    rare_colors <- names(color_counts)[color_counts < min_color_count]
    if (length(rare_colors) > 0) {
      node_df$.color[node_df$.color %in% rare_colors] <- "Other"
      cat("Merged into 'Other' (phylum, n<", min_color_count, "):", paste(
          rare_colors, collapse = ", "), "\n")
    }

    node_df$.degree <- suppressWarnings(as.numeric(node_df[[size_attr]]))
    node_df$.degree[is.na(node_df$.degree)] <- 0

    node_df$.label <- if (is.null(label_rank)) node_df$name else as.character(
      node_df[[label_rank]])

    groups <- unique(node_df$.group)
    edges_hier <- dplyr::bind_rows(
      tibble(from = "origin", to = groups),
      tibble(from = node_df$.group, to = node_df$name)
    )

    vertices <- dplyr::bind_rows(
      tibble(name = "origin", .color = NA_character_, .group = NA_character_,
        .label = NA_character_, .degree = NA_real_),
      tibble(name = groups, .color = NA_character_, .group = groups,
        .label = groups, .degree = NA_real_),
      node_df %>% dplyr::transmute(name, .color, .group, .label, .degree)
    )

    mygraph <- igraph::graph_from_data_frame(edges_hier, vertices = vertices)
    layout <- create_layout(mygraph, layout = "dendrogram", circular = TRUE)

    name_to_idx <- setNames(seq_len(nrow(layout)), layout$name)
    orig_edges <- igraph::as_data_frame(ig, what = "edges") %>%
      dplyr::select(from, to, sign) %>%
      dplyr::mutate(from_idx = name_to_idx[from], to_idx = name_to_idx[to]) %>%
      dplyr::filter(!is.na(from_idx), !is.na(to_idx))

    cat("Bundled edges plotted:", nrow(orig_edges), "of", ecount(ig),
      "original network edges\n")

    grp_pal <- make_palette(node_df$.group, group_colors)
    col_pal <- make_palette(node_df$.color, node_colors)

    n_color_levels <- length(col_pal)
    if (n_color_levels > 25 && show_node_legend) {
      message(color_attr, " has ", n_color_levels,
        " levels — legend suppressed for readability.")
      show_node_legend <- FALSE
    }

    leaf_layout <- layout %>%
      dplyr::filter(leaf) %>%
      dplyr::mutate(theta = atan2(y, x))

    degree_lookup <- setNames(node_df$.degree, node_df$name)
    leaf_layout$.degree_for_labels <- degree_lookup[leaf_layout$name]
    group_lookup <- setNames(node_df$.group, node_df$name)
    leaf_layout$.group <- group_lookup[leaf_layout$name]

    angle_raw <- leaf_layout$theta * 180 / pi
    leaf_layout$hjust <- ifelse(angle_raw < -90 | angle_raw > 90, 1, 0)
    leaf_layout$label_angle <- ifelse(angle_raw < -90 | angle_raw > 90,
      angle_raw + 180, angle_raw)

    labelable <- if (exclude_unclassified_labels) {
      leaf_layout %>% dplyr::filter(!grepl("^Unclassified", .label))
    } else {
      leaf_layout
    }

    top_leaf_names <- switch(
      label_mode,
      "all" = labelable$name,
      "top_n" = labelable %>% dplyr::arrange(desc(.degree_for_labels)) %>%
        dplyr::slice_head(n = label_top_n) %>% dplyr::pull(name),
      "per_group_top_n" = labelable %>% dplyr::group_by(.group) %>%
        dplyr::arrange(desc(.degree_for_labels), .by_group = TRUE) %>%
        dplyr::slice_head(n = label_top_n) %>% dplyr::ungroup() %>% dplyr::pull(
        name),
      stop("label_mode must be 'all', 'top_n', or 'per_group_top_n'")
    )

    label_df <- leaf_layout %>% dplyr::filter(name %in% top_leaf_names)
    cat("Labeling", nrow(label_df), "of", nrow(leaf_layout),
      "leaf nodes (mode:", label_mode, ")\n")

    # multiple, non-adjacent clusters around the ring — each gets its OWN
    # band polygon (correct boundaries), while sharing one fill color and
    # one legend entry via .group.
    leaf_sorted <- leaf_layout %>% dplyr::arrange(theta)
    n_leaf <- nrow(leaf_sorted)

    theta_vec <- leaf_sorted$theta
    next_theta <- c(theta_vec[-1], theta_vec[1] + 2 * pi)
    gap_to_next <- next_theta - theta_vec
    leaf_sorted$right_bound <- theta_vec + gap_to_next / 2

    prev_theta <- c(theta_vec[n_leaf] - 2 * pi, theta_vec[-n_leaf])
    gap_to_prev <- theta_vec - prev_theta
    leaf_sorted$left_bound <- theta_vec - gap_to_prev / 2

    run_lengths <- rle(leaf_sorted$.group)$lengths
    leaf_sorted$run_id <- rep(seq_along(run_lengths), times = run_lengths)

    n_runs <- max(leaf_sorted$run_id)
    if (n_runs > 1 && leaf_sorted$.group[leaf_sorted$run_id == 1][1] ==
      leaf_sorted$.group[leaf_sorted$run_id == n_runs][1]) {
      leaf_sorted$run_id[leaf_sorted$run_id == n_runs] <- 1
    }

    run_bounds <- leaf_sorted %>%
      dplyr::group_by(run_id) %>%
      dplyr::summarise(
        .group = dplyr::first(.group),
        theta_min = dplyr::first(left_bound),
        theta_max = dplyr::last(right_bound),
        .groups = "drop"
      ) %>%
      dplyr::filter(!is.na(.group))

    build_band_polygon <- function(theta_min, theta_max, group_name, run_id,
      n_points = 60) {
      thetas <- seq(theta_min, theta_max, length.out = n_points)
      outer_pts <- tibble(x = band_outer_r * cos(thetas), y = band_outer_r *
        sin(thetas))
      inner_pts <- tibble(x = band_inner_r * cos(rev(thetas)), y =
          band_inner_r * sin(rev(thetas)))
      tibble(x = c(outer_pts$x, inner_pts$x), y = c(outer_pts$y, inner_pts$y),
        .group = group_name, run_id = run_id)
    }

    band_polygons <- purrr::pmap_dfr(
      list(run_bounds$theta_min, run_bounds$theta_max, run_bounds$.group,
        run_bounds$run_id),
      build_band_polygon
    )

    cat("Guild bands drawn as", nrow(run_bounds), "contiguous segments across",
      length(unique(run_bounds$.group)), "guild categories\n")

    p <- ggraph(layout) +
      geom_polygon(
        data = band_polygons, aes(x = x, y = y, fill = .group, group = run_id),
        alpha = band_alpha, color = "white", linewidth = 0.3
      ) +
      scale_fill_manual(values = grp_pal, name = group_attr, guide =
          guide_legend(order = 1)) +
      new_scale_fill() +
      geom_conn_bundle(
        data = get_con(from = orig_edges$from_idx, to = orig_edges$to_idx,
          sign = orig_edges$sign),
        aes(edge_colour = sign), alpha = edge_alpha, width = edge_width,
        tension = tension
      ) +
      scale_edge_colour_manual(values = c(positive = "blue", negative = "red"),
        name = "Association sign") +
      ggnewscale::new_scale_color() +
      geom_node_point(
        data = function(d) dplyr::filter(d, leaf),
        aes(color = .color), size = node_size, alpha = 1
      ) +
      scale_color_manual(values = col_pal, name = color_attr,
        guide = if (show_node_legend) guide_legend(order = 3) else "none") +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(x = x * 1.05, y = y * 1.05, label = .label),
        size = label_size, fontface = "italic", color = "black",
        segment.size = 0.15, segment.color = "black",
        max.overlaps = Inf, box.padding = 0.15, point.padding = 0.05,
        min.segment.length = 0, seed = 123
      ) +
      coord_fixed(clip = "off") +
      expand_limits(x = c(-1.8, 1.8), y = c(-1.8, 1.8)) +
      labs(title = veg_group) +
      theme_void(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        legend.text = element_text(size = 7.5), legend.title = element_text(
          size = 8.5))

    p
  }

  # Rebuild all four panels with LARGER internal elements — bigger nodes,
  # bigger labels, thicker band/edge lines — so they don't look small once
  # squeezed into a 2x2 grid
  Plot_SMF <- plot_circular_bundle(
    veg_networks$SMF, "SMF",
    group_attr = "Guild_class", color_attr = "Phylum", size_attr = "degree",
    label_mode = "per_group_top_n", label_top_n = 10,
    # was 3.4 — bigger text
    label_size = 5.5,
    # NEW — bigger node points (if your function exposes this)
    node_size = 2.5,
    # thicker bundled edges
    edge_width = 0.4,
    band_alpha = 0.9,
    group_colors = guild_seed_colors,
    node_colors  = phylum_seed_colors
  )

  Plot_MF <- plot_circular_bundle(
    veg_networks$MF, "MF",
    group_attr = "Guild_class", color_attr = "Phylum", size_attr = "degree",
    label_mode = "per_group_top_n", label_top_n = 10,
    label_size = 5.5,
    node_size = 2.5,
    edge_width = 0.4,
    band_alpha = 0.9,
    group_colors = guild_seed_colors,
    node_colors  = phylum_seed_colors
  )

  Plot_UMF <- plot_circular_bundle(
    veg_networks$UMF, "UMF",
    group_attr = "Guild_class", color_attr = "Phylum", size_attr = "degree",
    label_mode = "per_group_top_n", label_top_n = 10,
    label_size = 5.5,
    node_size = 2.5,
    edge_width = 0.4,
    band_alpha = 0.9,
    group_colors = guild_seed_colors,
    node_colors  = phylum_seed_colors
  )

  Plot_CR <- plot_circular_bundle(
    veg_networks$CR, "CR",
    group_attr = "Guild_class", color_attr = "Phylum", size_attr = "degree",
    label_mode = "per_group_top_n", label_top_n = 10,
    label_size = 5.5,
    node_size = 2.5,
    edge_width = 0.4,
    band_alpha = 0.9,
    group_colors = guild_seed_colors,
    node_colors  = phylum_seed_colors
  )

  # Strip title/legend, combine as before
  Plot_SMF_clean <- Plot_SMF + labs(title = NULL) + theme(legend.position =
      "none")
  Plot_MF_clean  <- Plot_MF  + labs(title = NULL) + theme(legend.position =
      "none")
  Plot_UMF_clean <- Plot_UMF + labs(title = NULL) + theme(legend.position =
      "none")
  Plot_CR_clean  <- Plot_CR  + labs(title = NULL) + theme(legend.position =
      "none")

  combined_circular_bundles <- (Plot_SMF_clean | Plot_MF_clean) /
    (Plot_UMF_clean | Plot_CR_clean) +
  # force equal panel sizes
    plot_layout(widths = c(1, 1), heights = c(1, 1))

  combined_circular_bundles

  ggsave("circularbundle_combined_ALL_2.png", combined_circular_bundles,
    width = 40, height = 38, dpi = 300, bg = "white", limitsize = FALSE)
  ggsave("circularbundle_combined_ALL_2.pdf", combined_circular_bundles,
    device = cairo_pdf,
    width = 40, height = 38, bg = "white")
  # Create a plot satandlone legend
  # Check whether all four networks share the same categories first
  all_guilds <- unique(c(
      V(veg_networks$SMF$graph)$Guild_class, V(
        veg_networks$MF$graph)$Guild_class,
      V(veg_networks$UMF$graph)$Guild_class, V(
        veg_networks$CR$graph)$Guild_class
  ))
  all_phyla <- unique(c(
      V(veg_networks$SMF$graph)$Phylum, V(veg_networks$MF$graph)$Phylum,
      V(veg_networks$UMF$graph)$Phylum, V(veg_networks$CR$graph)$Phylum
  ))

  cat("Guild categories across ALL four groups:", paste(sort(all_guilds),
      collapse = ", "), "\n")
  cat("Phylum categories across ALL four groups:", paste(sort(all_phyla),
      collapse = ", "), "\n")

  cat("\nGuild categories present in SMF's legend:",
    paste(sort(unique(V(veg_networks$SMF$graph)$Guild_class)), collapse = ", "),
    "\n")

  # Build a minimal dummy dataset containing EVERY guild and phylum category
  # across all four groups — used only to generate a complete legend, not
  # to plot any real data points
  all_guilds <- c("Animal pathogen", "Dung saprotroph", "Endophyte",
    "Fungal parasite",
    "Generalist pathogen", "Generalist saprotroph", "Lichenized", "Mixed",
    "Mycorrhizal", "Nematophagous", "Plant pathogen", "Plant saprotroph",
    "Unassigned", "Wood saprotroph")

  all_phyla <- c("Aphelidiomycota", "Ascomycota", "Basidiobolomycota",
    "Basidiomycota",
    "Chytridiomycota", "Glomeromycota", "Kickxellomycota", "Monoblepharomycota",
    "Mortierellomycota", "Mucoromycota", "Rozellomycota", "Unclassified",
    "Zoopagomycota")

  dummy_guild_df <- tibble(x = seq_along(all_guilds), y = 1, Guild_class =
      all_guilds)
  dummy_phylum_df <- tibble(x = seq_along(all_phyla), y = 1, Phylum = all_phyla)
  dummy_sign_df <- tibble(x = 1:2, y = 1, sign = c("positive", "negative"))

  # Build three small invisible dummy plots, one per legend, using your
  # EXACT seed palettes — invisible geoms (size = 0) since we only want
  # the legend keys, not visible points
  p_dummy_guild <- ggplot(dummy_guild_df, aes(x = x, y = y, fill =
        Guild_class)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = guild_seed_colors, name = "Guild_class") +
    theme_void()

  p_dummy_phylum <- ggplot(dummy_phylum_df, aes(x, y, color = Phylum)) +
    geom_point(size = 0) +
    scale_color_manual(values = phylum_seed_colors, name = "Phylum",
      guide = guide_legend(override.aes = list(size = 4))) +
    theme_void()

  p_dummy_sign <- ggplot(dummy_sign_df, aes(x, y, color = sign)) +
    geom_line(aes(group = sign), linewidth = 0) +
    scale_color_manual(values = c(positive = "blue", negative = "red"),
      name = "Association sign",
      guide = guide_legend(override.aes = list(linewidth = 2))) +
    theme_void()

  legend_guild  <- get_legend(p_dummy_guild)
  legend_phylum <- get_legend(p_dummy_phylum)
  legend_sign   <- get_legend(p_dummy_sign)

  combined_legend <- plot_grid(
    legend_guild, legend_phylum, legend_sign,
    ncol = 1, align = "v",
    rel_heights = c(1.4, 1.3, 0.5)
  )

  combined_legend

  ggsave("network_legend_complete.png", combined_legend, width = 4, height = 12,
    dpi = 300, bg = "white")
  save_pdf("network_legend_complete.pdf", combined_legend, width = 4, height =
      12)
  # Export CSV tables + network statistics for verification

  export_network_results <- function(net_obj, veg_group,
    group_attr = "Guild_class",
    color_attr = "Phylum",
    min_group_count = 1,
    min_color_count = 1) {

    ig <- net_obj$graph
    node_df <- igraph::as_data_frame(ig, what = "vertices")

    # 1. Node-level table — RAW attribute values (before any "Other"
    # merging), so you can directly verify every genus's true
    # Guild_class/Phylum/degree as it exists in the network
    node_table <- node_df %>%
      dplyr::transmute(
        Genus = name,
        Guild_class_raw = .data[[group_attr]],
        Phylum_raw = .data[[color_attr]],
        Degree = degree(ig)[match(name, V(ig)$name)]
      ) %>%
      dplyr::arrange(desc(Degree))

    # 2. Apply the SAME merge-to-"Other" logic the plot uses, as an
    # additional column — lets you directly compare raw vs. what
    # actually got plotted/colored
    group_counts <- table(node_table$Guild_class_raw)
    rare_groups <- names(group_counts)[group_counts < min_group_count]
    node_table$Guild_class_plotted <- ifelse(node_table$Guild_class_raw %in%
        rare_groups,
      "Other", node_table$Guild_class_raw)

    color_counts <- table(node_table$Phylum_raw)
    rare_colors <- names(color_counts)[color_counts < min_color_count]
    node_table$Phylum_plotted <- ifelse(node_table$Phylum_raw %in% rare_colors,
      "Other", node_table$Phylum_raw)

    write.csv(node_table, paste0("network_nodes_", veg_group, ".csv"),
      row.names = FALSE)
    cat("Exported: network_nodes_", veg_group, ".csv (", nrow(node_table),
      " genera)\n", sep = "")

    # 3. Edge-level table — every co-occurrence relationship, with sign
    # and weight, plus both endpoints' guild/phylum for easy filtering
    edge_table <- igraph::as_data_frame(ig, what = "edges") %>%
      dplyr::left_join(node_table %>% dplyr::select(Genus, Guild_class_raw,
          Phylum_raw),
        by = c("from" = "Genus")) %>%
      dplyr::rename(from_Guild = Guild_class_raw, from_Phylum = Phylum_raw) %>%
      dplyr::left_join(node_table %>% dplyr::select(Genus, Guild_class_raw,
          Phylum_raw),
        by = c("to" = "Genus")) %>%
      dplyr::rename(to_Guild = Guild_class_raw, to_Phylum = Phylum_raw)

    write.csv(edge_table, paste0("network_edges_", veg_group, ".csv"),
      row.names = FALSE)
    cat("Exported: network_edges_", veg_group, ".csv (", nrow(edge_table),
      " edges)\n", sep = "")

    # 4. Guild/Phylum summary counts — the exact numbers behind the bands
    # and legend, so you can directly cross-check what SHOULD appear
    guild_summary <- node_table %>% dplyr::count(Guild_class_raw,
      Guild_class_plotted, sort = TRUE)
    phylum_summary <- node_table %>% dplyr::count(Phylum_raw, Phylum_plotted,
      sort = TRUE)

    write.csv(guild_summary, paste0("guild_summary_", veg_group, ".csv"),
      row.names = FALSE)
    write.csv(phylum_summary, paste0("phylum_summary_", veg_group, ".csv"),
      row.names = FALSE)

    # 5. Overall network topology statistics
    clustering <- cluster_fast_greedy(ig, weights = abs(E(ig)$weight))

    topology_stats <- tibble(
      Veg = veg_group,
      n_nodes = vcount(ig),
      n_edges = ecount(ig),
      n_positive_edges = sum(E(ig)$sign == "positive"),
      n_negative_edges = sum(E(ig)$sign == "negative"),
      avg_degree = round(mean(degree(ig)), 2),
      edge_density = round(edge_density(ig), 4),
      transitivity = round(transitivity(ig, type = "global"), 4),
      modularity = round(modularity(clustering), 4),
      n_modules = length(unique(membership(clustering))),
      stars_stability = round(net_obj$stability, 4),
      n_guild_categories = length(unique(node_table$Guild_class_raw)),
      n_phylum_categories = length(unique(node_table$Phylum_raw))
    )

    write.csv(topology_stats, paste0("network_topology_stats_", veg_group,
        ".csv"), row.names = FALSE)

    cat("\n=== Topology stats for", veg_group, "===\n")
    print(topology_stats)

    cat("\n=== Guild summary (raw vs. plotted) for", veg_group, "===\n")
    print(guild_summary)

    cat("\n=== Phylum summary (raw vs. plotted) for", veg_group, "===\n")
    print(phylum_summary)

    invisible(list(nodes = node_table, edges = edge_table,
        guild_summary = guild_summary, phylum_summary = phylum_summary,
        topology = topology_stats))
  }

  # Run for SMF
  smf_results <- export_network_results(veg_networks$SMF, "SMF")
  mf_results  <- export_network_results(veg_networks$MF,  "MF")
  umf_results <- export_network_results(veg_networks$UMF, "UMF")

  # ----------------------------------------------------------------------------
  # 12. Network robustness simulations
  # ----------------------------------------------------------------------------
  cr_results  <- export_network_results(veg_networks$CR,  "CR")
  # Network robustness/vulnerability metrics — targeted vs. random node removal

  simulate_robustness <- function(ig, n_sims = 50) {
    n <- vcount(ig)

    # --- Targeted attack: remove highest-degree node each step ---
    ig_target <- ig
    target_curve <- numeric(n)
    for (i in seq_len(n)) {
      target_curve[i] <- max(components(ig_target)$csize) / n
      if (vcount(ig_target) == 0) break
      deg <- degree(ig_target)
      remove_node <- names(which.max(deg))
      ig_target <- delete_vertices(ig_target, remove_node)
    }

    # --- Random attack: averaged over n_sims random removal orders ---
    random_curves <- matrix(NA, nrow = n, ncol = n_sims)
    for (s in seq_len(n_sims)) {
      ig_rand <- ig
      order_remove <- sample(V(ig)$name)
      for (i in seq_len(n)) {
        random_curves[i, s] <- max(components(ig_rand)$csize) / n
        if (vcount(ig_rand) == 0) break
        ig_rand <- delete_vertices(ig_rand, order_remove[i])
      }
    }
    random_curve <- rowMeans(random_curves, na.rm = TRUE)

    tibble(prop_removed = seq_len(n) / n, targeted = target_curve, random =
        random_curve)
  }

  # Summarize each curve into reportable single-number metrics:
  # - AUC (area under the fragmentation curve): higher = more robust,
  # since the network stays largely connected across more removals
  # - Critical fraction: proportion of nodes removed at which the
  # largest connected component first drops below 50% of the network
  # - Vulnerability index: AUC(random) - AUC(targeted), a common summary
  # of how much MORE damage targeted (hub-based) attack does relative
  # to random loss — larger gap = more hub-dependent/fragile structure
  summarize_robustness <- function(robustness_df, veg_group) {
    auc <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

    auc_targeted <- auc(robustness_df$prop_removed, robustness_df$targeted)
    auc_random   <- auc(robustness_df$prop_removed, robustness_df$random)

    critical_targeted <- robustness_df$prop_removed[which(
        robustness_df$targeted < 0.5)[1]]
    critical_random   <- robustness_df$prop_removed[which(robustness_df$random <
        0.5)[1]]

    tibble(
      Veg = veg_group,
      AUC_targeted = round(auc_targeted, 4),
      AUC_random = round(auc_random, 4),
      Vulnerability_index = round(auc_random - auc_targeted, 4),
      Critical_fraction_targeted = round(critical_targeted, 4),
      Critical_fraction_random = round(critical_random, 4)
    )
  }

  # Run for SMF, export both the full curve data and the summary metrics
  set.seed(123)
  robustness_curve_SMF <- simulate_robustness(veg_networks$SMF$graph)
  write.csv(robustness_curve_SMF, "robustness_curve_SMF.csv", row.names = FALSE)

  robustness_summary_SMF <- summarize_robustness(robustness_curve_SMF, "SMF")
  write.csv(robustness_summary_SMF, "robustness_summary_SMF.csv", row.names =
      FALSE)

  cat("=== Robustness summary for SMF ===\n")
  print(robustness_summary_SMF)

  set.seed(123)
  robustness_curve_MF <- simulate_robustness(veg_networks$MF$graph)
  robustness_summary_MF <- summarize_robustness(robustness_curve_MF, "MF")

  set.seed(123)
  robustness_curve_UMF <- simulate_robustness(veg_networks$UMF$graph)
  robustness_summary_UMF <- summarize_robustness(robustness_curve_UMF, "UMF")

  set.seed(123)
  robustness_curve_CR <- simulate_robustness(veg_networks$CR$graph)
  robustness_summary_CR <- summarize_robustness(robustness_curve_CR, "CR")

  write.csv(robustness_curve_MF, "robustness_curve_MF.csv", row.names = FALSE)
  write.csv(robustness_curve_UMF, "robustness_curve_UMF.csv", row.names = FALSE)
  write.csv(robustness_curve_CR, "robustness_curve_CR.csv", row.names = FALSE)

  write.csv(robustness_summary_MF, "robustness_summary_MF.csv", row.names =
      FALSE)
  write.csv(robustness_summary_UMF, "robustness_summary_UMF.csv", row.names =
      FALSE)
  write.csv(robustness_summary_CR, "robustness_summary_CR.csv", row.names =
      FALSE)

  # Combined summary table across all four groups — the key comparison table
  robustness_all <- dplyr::bind_rows(robustness_summary_SMF,
    robustness_summary_MF,
    robustness_summary_UMF, robustness_summary_CR)
  write.csv(robustness_all, "robustness_summary_ALL_GROUPS.csv", row.names =
      FALSE)

  cat("\n=== Robustness summary across all vegetation types ===\n")
  print(robustness_all)
  # Plot 1: Fragmentation curves — targeted vs random, all four groups overlaid

  veg_colors <- c(SMF = "#3B75AF", MF = "#519E3E", UMF = "#EF8A17", CR =
      "#C53A32")

  robustness_curves_all <- dplyr::bind_rows(
    robustness_curve_SMF %>% dplyr::mutate(Veg = "SMF"),
    robustness_curve_MF  %>% dplyr::mutate(Veg = "MF"),
    robustness_curve_UMF %>% dplyr::mutate(Veg = "UMF"),
    robustness_curve_CR  %>% dplyr::mutate(Veg = "CR")
  )

  robustness_long <- robustness_curves_all %>%
    pivot_longer(cols = c(targeted, random), names_to = "Attack_type",
      values_to = "Largest_component_frac") %>%
    dplyr::mutate(Veg = factor(Veg, levels = names(veg_colors)),
      Attack_type = factor(Attack_type, levels = c("targeted", "random"),
        labels = c("Targeted (highest-degree first)", "Random")))

  p_curves <- ggplot(robustness_long, aes(x = prop_removed, y =
        Largest_component_frac, color = Veg)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50",
      linewidth = 0.4) +
    facet_wrap(~ Attack_type) +
    scale_color_manual(values = veg_colors) +
    labs(
      x = "Proportion of nodes removed",
      y = "Largest connected component (proportion of network)",
      color = "Vegetation type",
      title = "Network fragmentation under targeted vs. random node removal",
      subtitle = "Dashed line = 50% connectivity threshold (critical fraction)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      legend.position = "right"
    )

  p_curves
  ggsave("robustness_fragmentation_curves.png", p_curves, width = 12, height =
      6, dpi = 300, bg = "white")
  ggsave("robustness_fragmentation_curves.pdf", p_curves, width = 12, height =
      6, bg = "white")

  # Plot 2: Summary metrics — grouped bar chart across all four groups
  robustness_summary_long <- robustness_all %>%
    dplyr::select(Veg, AUC_targeted, AUC_random, Critical_fraction_targeted,
      Critical_fraction_random) %>%
    pivot_longer(cols = -Veg, names_to = "Metric", values_to = "Value") %>%
    dplyr::mutate(
      Metric_type = ifelse(grepl("^AUC", Metric), "Area Under Curve (AUC)",
        "Critical Fraction"),
      Attack_type = ifelse(grepl("targeted", Metric), "Targeted", "Random"),
      Veg = factor(Veg, levels = names(veg_colors))
    )

  p_summary_bars <- ggplot(robustness_summary_long, aes(x = Veg, y = Value,
      fill = Attack_type)) +
    geom_col(position = "dodge", width = 0.7, color = "white", linewidth =
        0.2) +
    facet_wrap(~ Metric_type, scales = "free_y") +
    scale_fill_manual(values = c(Targeted = "red", Random = "blue")) +
    labs(x = "Vegetation type", y = "Value", fill = "Attack type",
      title = "Robustness summary metrics across vegetation types") +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 14)
    )

  p_summary_bars
  ggsave("robustness_summary_bars.png", p_summary_bars, width = 11, height = 6,
    dpi = 300, bg = "white")
  ggsave("robustness_summary_bars.pdf", p_summary_bars, width = 11, height = 6,
    bg = "white")

  # Plot 3: Vulnerability index alone — the single most important comparison
  p_vulnerability <- ggplot(robustness_all, aes(x = reorder(Veg, -
        Vulnerability_index), y = Vulnerability_index, fill = Veg)) +
    geom_col(width = 0.6, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = veg_colors) +
    labs(x = "Vegetation type", y =
        "Vulnerability index (AUC random \u2212 AUC targeted)",
      title = "Network vulnerability to targeted hub removal",
      subtitle = paste0("Higher values indicate greater reliance ",
        "on hub taxa for network integrity")) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9.5, color = "grey40")
    )

  p_vulnerability
  ggsave("robustness_vulnerability_index.png", p_vulnerability, width = 8,
    height = 6, dpi = 300, bg = "white")
  ggsave("robustness_vulnerability_index.pdf", p_vulnerability, width = 8,
    height = 6, bg = "white")
  # Optional: combine all three into one composite figure
  # Strip individual titles/subtitles from each panel before combining —
  # the combined figure's plot_annotation() already carries the overall title
  p_curves_clean <- p_curves + labs(title = NULL, subtitle = NULL)
  p_summary_bars_clean <- p_summary_bars + labs(title = NULL, subtitle = NULL)
  p_vulnerability_clean <- p_vulnerability + labs(title = NULL, subtitle = NULL)

  combined_robustness <- (p_curves_clean) / (p_summary_bars_clean |
      p_vulnerability_clean) +
    plot_annotation(tag_levels = "a", tag_suffix = ")")

  combined_robustness
  ggsave("robustness_combined_figure.png", combined_robustness, width = 15,
    height = 12, dpi = 300, bg = "white")
  ggsave("robustness_combined_figure.pdf", combined_robustness, device =
      cairo_pdf, width = 15, height = 12, bg = "white")

  # ----------------------------------------------------------------------------
  # 13. Within-module and among-module node roles
  # ----------------------------------------------------------------------------
  # Keystone taxa identification via Zi-Pi (within-module degree /
  # among-module connectivity) classification

  compute_zi_pi <- function(ig) {
    # Modularity-based community detection (weights must be non-negative,
    # so use absolute association strength — same convention used throughout)
    modules <- cluster_fast_greedy(ig, weights = abs(E(ig)$weight))
    mem <- membership(modules)

    zi <- numeric(vcount(ig))
    pi_val <- numeric(vcount(ig))
    names(zi) <- names(pi_val) <- V(ig)$name

    # --- Zi: within-module degree z-score ---
    for (m in unique(mem)) {
      module_nodes <- names(mem)[mem == m]
      sub_deg <- sapply(module_nodes, function(n) {
          neighbors_n <- neighbors(ig, n)
          sum(names(neighbors_n) %in% module_nodes)
      })
      mean_k <- mean(sub_deg)
      sd_k <- sd(sub_deg)
      zi[module_nodes] <- if (is.finite(sd_k) && sd_k > 0) (sub_deg - mean_k) /
        sd_k else 0
    }

    # --- Pi: among-module connectivity (participation coefficient) ---
    for (v in V(ig)$name) {
      neighbors_v <- neighbors(ig, v)
      k_total <- length(neighbors_v)
      if (k_total == 0) { pi_val[v] <- 0; next }
      module_counts <- table(mem[names(neighbors_v)])
      pi_val[v] <- 1 - sum((module_counts / k_total)^2)
    }

    tibble(
      Taxon = V(ig)$name,
      Zi = zi[V(ig)$name],
      Pi = pi_val[V(ig)$name],
      Module = mem[V(ig)$name],
      Degree = degree(ig)[V(ig)$name],
      Guild_class = V(ig)$Guild_class,
      Phylum = V(ig)$Phylum
    ) %>%
      dplyr::mutate(
      # Standard Guimerà-Amaral/Olesen thresholds — the accepted
      # convention in the network ecology literature (not arbitrary)
        Node_role = case_when(
          Zi > 2.5 & Pi > 0.62 ~ "Network hub",
          Zi > 2.5 & Pi <= 0.62 ~ "Module hub",
          Zi <= 2.5 & Pi > 0.62 ~ "Connector",
          TRUE ~ "Peripheral"
        ),
      # the taxa structurally important enough to matter for network
      # integrity, per the Zi-Pi framework
        Is_keystone = Node_role != "Peripheral"
      )
  }

  # Run across all four vegetation-type networks
  zi_pi_all <- purrr::imap_dfr(veg_networks, function(net_obj, veg_group) {
      compute_zi_pi(net_obj$graph) %>% dplyr::mutate(Veg = veg_group)
  })

  write.csv(zi_pi_all, "zi_pi_all_taxa.csv", row.names = FALSE)

  keystone_taxa <- zi_pi_all %>%
    dplyr::filter(Is_keystone) %>%
    dplyr::arrange(Veg, desc(Zi), desc(Pi))

  write.csv(keystone_taxa, "keystone_taxa_all_groups.csv", row.names = FALSE)

  cat("=== Keystone taxa counts by vegetation type and role ===\n")
  print(zi_pi_all %>% dplyr::count(Veg, Node_role))

  cat("\n=== Keystone taxa (all roles except Peripheral) ===\n")
  print(keystone_taxa %>% dplyr::select(Veg, Taxon, Zi, Pi, Node_role,
      Guild_class, Phylum), n = 30)

  # Guild composition of keystone taxa — are hubs disproportionately one
  # functional type per vegetation group? (the biologically interesting part)
  cat(paste0("\n=== Functional guild composition of ",
      "keystone taxa, by vegetation type ===\n"))
  print(keystone_taxa %>% dplyr::count(Veg, Guild_class) %>% dplyr::arrange(Veg,
      desc(n)))

  # Zi-Pi scatter plot — the standard visualization for this framework
  p_zipi <- ggplot(zi_pi_all, aes(x = Pi, y = Zi)) +
    geom_point(
      aes(fill = Node_role),
      shape = 21,
      colour = "black",
      size = 4,
      stroke = 0.25,
      alpha = 0.7
    )+
    geom_hline(yintercept = 2.5, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0.62, linetype = "dashed", color = "grey40") +
    facet_wrap(~ Veg) +
    scale_fill_manual(values = c(
        "Network hub" = "#B2182B",
        "Module hub"  = "orange",
        "Connector"   = "blue",
        "Peripheral"  = "green"
    ))+
    labs(x = "Among-module connectivity (Pi)", y = "Within-module degree (Zi)",
      fill = "Node role")+
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold"), plot.title = element_text(
        face = "bold"))

  p_zipi
  ggsave("zi_pi_classification_plot.png", p_zipi, width = 11, height = 9, dpi =
      300, bg = "white")
  ggsave("zi_pi_classification_plot.pdf", p_zipi, device = cairo_pdf, width =
      11, height = 9, bg = "white")

  cat("=== Exact Node_role counts per vegetation type ===\n")
  print(zi_pi_all %>% dplyr::count(Veg, Node_role) %>% dplyr::arrange(Veg,
      Node_role))

  cat(paste0("\n=== Which vegetation types have ZERO ",
      "module hubs or network hubs? ===\n"))
  zi_pi_all %>%
    dplyr::count(Veg, Node_role) %>%
    tidyr::pivot_wider(names_from = Node_role, values_from = n, values_fill =
        0) %>%
    print()

  # Keystone taxa counts by vegetation type and network role

  keystone_counts_summary <- zi_pi_all %>%
    dplyr::count(Veg, Node_role) %>%
    tidyr::pivot_wider(names_from = Node_role, values_from = n, values_fill =
        0) %>%
  # Ensure consistent column order regardless of which roles happened to
  # appear first in the data
    dplyr::select(Veg, any_of(c("Network hub", "Module hub", "Connector",
          "Peripheral")))

  cat("=== Keystone taxa counts by vegetation type and network role ===\n")
  print(keystone_counts_summary)

  write.csv(keystone_counts_summary, "Table_keystone_counts_by_role.csv",
    row.names = FALSE)

  # ---------------------------------------------------------------------------
  # 14. Mantel tests and environmental correlations
  # ----------------------------------------------------------------------------
  # Combined correlation matrix + Mantel test figure (linkET-style)
  # Using ALL 17 environmental variables — no VIF screening needed for Mantel,
  # since each Mantel test is univariate/independent, unlike RDA's joint model
  # remotes::install_github("Hy4m/linkET", force = TRUE)   # run once if not
  # installed

  set.seed(123)

  # 1. Environmental variable matrix — ALL 17 candidate variables
  # (env_matrix, not env_reduced — Mantel doesn't require VIF screening)
  env_for_plot <- env_matrix

  # 2. Build distance matrices for each "response" metric
  hill_matched <- hill_df %>%
    dplyr::filter(Sample %in% rownames(env_for_plot)) %>%
    dplyr::arrange(match(Sample, rownames(env_for_plot)))

  q0_dist <- dist(hill_matched$q0, method = "euclidean")
  q1_dist <- dist(hill_matched$q1, method = "euclidean")
  q2_dist <- dist(hill_matched$q2, method = "euclidean")

  community_dist <- as.dist(as.matrix(aitchison_dist)[rownames(env_for_plot),
      rownames(env_for_plot)])

  response_list <- list(
    Richness_q0 = q0_dist,
    Shannon_q1  = q1_dist,
    Simpson_q2  = q2_dist,
    Community   = community_dist
  )

  # 3. Build a distance matrix for EACH of the 17 environmental variables
  env_dist_list <- map(colnames(env_for_plot), function(varname) {
      dist(scale(env_for_plot[[varname]]), method = "euclidean")
  })
  names(env_dist_list) <- colnames(env_for_plot)

  # 4. Loop every response distance matrix against every one of the 17
  # environmental variables' distance matrices
  mantel_results <- map_dfr(names(response_list), function(resp_name) {
      map_dfr(names(env_dist_list), function(env_name) {
          test <- vegan::mantel(response_list[[resp_name]], env_dist_list[[
                env_name]],
            method = "spearman", permutations = 999)
          tibble(spec = resp_name, env = env_name, r = test$statistic, p =
              test$signif)
      })
  })

  mantel_results <- mantel_results %>%
    dplyr::mutate(
      r_category = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), labels = c("< 0.2",
          "0.2 - 0.4", ">= 0.4")),
      p_category = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), labels = c(
          "< 0.01", "0.01 - 0.05", ">= 0.05"))
    )

  cat(paste0("=== Mantel test results: ALL 17 ",
      "environmental variables, manual loop ===\n"))
  print(mantel_results, n = Inf)
  write.csv(mantel_results, "mantel_test_results_all_17_variables.csv",
    row.names = FALSE)

  # 5. Build the combined figure — 17-variable correlation matrix + Mantel
  # lines, matching the linkET reference figure convention
  # Add significance asterisks to the environmental-variable correlation
  # matrix squares (Pearson correlation, from qcorrplot)
  env_corr <- correlate(env_for_plot, method = "pearson")

  # 2. Extract p-values into a long-format table with significance stars
  p_mat <- env_corr$p
  r_mat <- env_corr$r

  sig_df <- as.data.frame(as.table(p_mat)) %>%
    dplyr::rename(Var1 = Var1, Var2 = Var2, p = Freq) %>%
    dplyr::left_join(
      as.data.frame(as.table(r_mat)) %>% dplyr::rename(Var1 = Var1, Var2 =
          Var2, r = Freq),
      by = c("Var1", "Var2")
    ) %>%
  # skip diagonal
    dplyr::filter(as.character(Var1) != as.character(Var2)) %>%
    dplyr::mutate(
      stars = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE ~ ""
      )
    )

  # 3. Rebuild the plot, adding a text layer for significance stars on top
  # of the correlation squares — geom_mark() is linkET's own text-overlay
  # function, designed specifically for this purpose
  # Correlation matrix squares with ONLY significance stars shown (no r
  # values)

  env_corr <- correlate(env_for_plot, method = "pearson")

  p_combined <- qcorrplot(env_corr, type = "lower", diag = FALSE) +
    geom_square() +
    geom_mark(only_mark = TRUE, sig_thres = 0.05, size = 3.2, colour =
        "black") +
    geom_couple(
      aes(colour = p_category, size = r_category),
      data = mantel_results,
      node = list(size = 3, color = "darkblue"),
      curvature = 0.15
    ) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#f90e06",
      name = "Pearson's r") +
    scale_size_manual(values = c(0.5, 1.2, 2.2), name = "Mantel's r") +
    scale_colour_manual(values = c("darkblue", "#00A651", "hotpink"), name =
        "Mantel's p")

  p_combined

  ggsave("mantel_correlation_combined_stars_only.pdf", p_combined, device =
      cairo_pdf, width = 12, height = 10, units = "in")
  ggsave("mantel_correlation_combined_stars_only.png", p_combined, width = 12,
    height = 10, dpi = 300, bg = "white")

  # TABLES accompanying the combined figures

  # Table A — Pairwise Pearson correlation matrix among all 17 environmental
  # variables
  env_corr_matrix <- cor(env_for_plot, method = "pearson", use =
      "pairwise.complete.obs")

  cat(paste0("\n=== Table A: Environmental variable ",
      "correlation matrix (Pearson's r, n=17) ", "===\n"))
  print(round(env_corr_matrix, 3))
  write.csv(round(env_corr_matrix, 3),
    "TableA_environmental_correlation_matrix_17vars.csv")

  # Table A (long format) — same correlations, tidy long table with
  # significance
  env_corr_test <- Hmisc::rcorr(as.matrix(env_for_plot), type = "pearson")

  env_corr_long <- as.data.frame(as.table(env_corr_test$r)) %>%
    dplyr::rename(Var1 = Var1, Var2 = Var2, r = Freq) %>%
    dplyr::left_join(
      as.data.frame(as.table(env_corr_test$P)) %>% dplyr::rename(Var1 = Var1,
        Var2 = Var2, p = Freq),
      by = c("Var1", "Var2")
    ) %>%
    dplyr::filter(Var1 != Var2) %>%
    dplyr::mutate(pair_id = paste(pmin(as.character(Var1), as.character(
            Var2)),
        pmax(as.character(Var1), as.character(Var2)))) %>%
    dplyr::distinct(pair_id, .keep_all = TRUE) %>%
    dplyr::select(Var1, Var2, r, p) %>%
    dplyr::arrange(desc(abs(r)))

  cat(paste0("\n=== Table A (long format): pairwise ",
      "correlations, sorted by |r|, n=17 ", "variables ===\n"))
  print(env_corr_long)
  write.csv(env_corr_long, "TableA_environmental_correlation_long_17vars.csv",
    row.names = FALSE)

  # Table B — Mantel test results (diversity metrics vs. ALL 17 environmental
  # variables), manuscript-ready formatting
  mantel_table_clean <- mantel_results %>%
    dplyr::transmute(
      `Diversity metric` = spec,
      `Environmental variable` = env,
      `Mantel's r` = round(r, 3),
      `p-value` = p,
      `Significance` = case_when(
        p < 0.01 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ "n.s."
      )
    ) %>%
    dplyr::arrange(`Diversity metric`, desc(`Mantel's r`))

  cat(paste0("\n=== Table B: Mantel test results, ALL ",
      "17 environmental variables ===\n"))
  print(mantel_table_clean, n = Inf)
  write.csv(mantel_table_clean, "TableB_mantel_test_results_17vars.csv",
    row.names = FALSE)

  # Table B (significant only) — compact main-text-ready version
  mantel_table_significant <- mantel_table_clean %>% dplyr::filter(`p-value` <
    0.05)

  cat("\n=== Table B (significant only, p<0.05) ===\n")
  print(mantel_table_significant, n = Inf)
  write.csv(mantel_table_significant,
    "TableB_mantel_significant_only_17vars.csv", row.names = FALSE)

  # Cross-reference: compare community-column Mantel results against your
  # earlier 17-variable envfit table, to directly verify the agreement/
  # disagreement pattern discussed (pH/CC/SpR strong in both; Cla/TN
  # borderline)
  mantel_vs_envfit <- mantel_table_clean %>%
    dplyr::filter(`Diversity metric` == "Community") %>%
    dplyr::select(`Environmental variable`, `Mantel's r`, `Mantel p` =
        `p-value`) %>%
    dplyr::left_join(
      envfit_pcoa_df %>% dplyr::transmute(`Environmental variable` = Variable,
        `envfit r2` = round(r2, 3), `envfit p` = p),
      by = "Environmental variable"
    ) %>%
    dplyr::arrange(desc(`envfit r2`))

  cat(paste0("\n=== Cross-reference: Mantel (community) ",
      "vs. envfit (PCoA), all 17 variables ===\n"))
  print(mantel_vs_envfit, n = Inf)
  write.csv(mantel_vs_envfit, "SuppTable_Mantel_vs_envfit_comparison.csv",
    row.names = FALSE)

  # ----------------------------------------------------------------------------
  # 15. Network topology summary
  # ----------------------------------------------------------------------------
  # Network topology summary — regenerated for the final genus-level
  # per-Veg-group networks (veg_networks)

  compute_topology <- function(net_obj, veg_group) {
    ig <- net_obj$graph

    # cluster_fast_greedy requires non-negative weights — use absolute
    # edge strength for community detection; sign stays a separate attribute
    clustering <- cluster_fast_greedy(ig, weights = abs(E(ig)$weight))

    tibble(
      `Vegetation type` = veg_group,
      Nodes = vcount(ig),
      Edges = ecount(ig),
      `Avg. degree` = round(mean(degree(ig)), 2),
      `Edge density` = round(edge_density(ig), 4),
      Modularity = round(modularity(clustering), 4),
      `Pos:Neg ratio` = round(sum(E(ig)$sign == "positive") / max(sum(E(
              ig)$sign == "negative"), 1), 2),
      `StARS stability` = round(net_obj$stability, 4)
    )
  }

  network_topology_summary <- imap_dfr(veg_networks, compute_topology)

  cat(paste0("=== Network topology summary ",
      "(genus-level, per vegetation/land-use ", "type) ===\n"))
  print(network_topology_summary)

  write.csv(network_topology_summary, "Table_network_topology_summary.csv",
    row.names = FALSE)

  # ----------------------------------------------------------------------------
  # 16. Supplementary taxonomic and guild tables
  # ----------------------------------------------------------------------------
  # Table (a): Resolved/Unresolved OTU counts + distinct resolved taxa per rank,
  # per vegetation/land-use type

  compute_resolution_table <- function(ps, rank) {
    tax_df <- as.data.frame(tax_table(ps))
    veg_vec <- as(sample_data(ps), "data.frame")$Veg
    otu_mat <- as(otu_table(ps), "matrix")
    if (!taxa_are_rows(ps)) otu_mat <- t(otu_mat)

    map_dfr(unique(veg_vec), function(v) {
        samples_v <- sample_names(ps)[veg_vec == v]
        present_otus <- rownames(otu_mat)[rowSums(otu_mat[, samples_v, drop =
                FALSE] > 0) > 0]
        tax_sub <- tax_df[present_otus, ]

        is_resolved <- !is.na(tax_sub[[rank]]) & tax_sub[[rank]] !=
        "Unclassified" & tax_sub[[rank]] != ""

        tibble(
          Veg = v, Rank = rank,
          Resolved_OTUs = sum(is_resolved),
          Unresolved_OTUs = sum(!is_resolved),
          Total_OTUs = length(present_otus),
          Distinct_resolved_taxa = length(unique(tax_sub[[rank]][is_resolved]))
        )
    })
  }

  resolution_table <- purrr::map_dfr(c("Phylum", "Class", "Order", "Family",
      "Genus", "Species"),
    ~ compute_resolution_table(fungi_physeq, .x))

  resolution_table <- resolution_table %>%
    dplyr::mutate(Rank = factor(Rank, levels = c("Phylum", "Class", "Order",
          "Family", "Genus", "Species"))) %>%
    dplyr::arrange(Veg, Rank)

  print(resolution_table, n = Inf)
  write.csv(resolution_table, "Table_taxonomic_resolution_summary.csv",
    row.names = FALSE)

  # Table (b): Mean relative abundance per functional guild, per vegetation/
  # land-use type
  fungi_rel <- transform_sample_counts(fungi_physeq, function(x) x / sum(x))
  df <- psmelt(fungi_rel)

  df$Guild_class[is.na(df$Guild_class) | df$Guild_class == ""] <- "Unassigned"

  guild_table <- df %>%
    dplyr::group_by(Sample, Veg, Guild_class) %>%
    dplyr::summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    dplyr::group_by(Veg, Guild_class) %>%
    dplyr::summarise(Mean_RA = mean(Abundance), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Veg, values_from = Mean_RA, values_fill =
        0) %>%
  # convert to %
    dplyr::mutate(across(where(is.numeric), ~ round(.x * 100, 2)))

  print(guild_table, n = Inf)
  write.csv(guild_table, "Table_guild_composition_by_Veg.csv", row.names =
      FALSE)

  writeLines(capture.output(sessionInfo()), "session_info.txt")
  message("Analysis complete. Outputs: ", normalizePath(getwd()))
  invisible(NULL)
}

run_fungi_analysis()
# ---------------------------------END OF SCRIPT--------------------------------
# //////////////////////////////////////////////////////////////////////////////
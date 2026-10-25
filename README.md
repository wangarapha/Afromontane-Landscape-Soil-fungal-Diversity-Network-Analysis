# Soil fungal diversity and association networks in the Bwindi landscape

This repository contains the R code used for the analyses, figures and supplementary tables in the manuscript:

>Authors: Rapheal Wangalwa, Jovan Ssemanda, Justine Nakintu1,3, Daisy Namirembe4, Desire Akanyijuka5,6, Tonny Chemutos1, Hilda Ikiriza1, Eunice A. Olet1 
Manuscript Title : When forest meets cropland: metabarcoding reveals shifts in soil fungal diversity and network structure across land use and elevation in an Afromontane landscape 
Submitted to *Landscape Ecology*.

The analysis works from a phyloseq object of soil fungal OTUs obtained by metabarcoding, with samples grouped into four habitat types (`SMF`, `MF`, `UMF` and `CR`; [HABITAT DEFINITIONS TO BE ADDED]).

## Data availability

The data will be made publicly available upon publication of the manuscript. The input file is therefore not included here. Once the manuscript is published, this section will be updated with the repository and accession number or DOI: [TO BE ADDED].

## Repository structure

```
.
├── R/
│   ├── Fungi_Paper_Script.R   # Full analysis pipeline
│   └── packages.R             # Installs the required R packages
├── data/
│   └── README.md              # Expected input format (data not included)
├── fungi-analysis.Rproj       # RStudio project; sets the working directory
├── .gitignore
└── README.md
```

## What the script does

`R/Fungi_Paper_Script.R` runs as a single function, `run_fungi_analysis()`, organised in numbered sections:

1. Configuration and input checks
2. Sequencing overview (reads per sample and per habitat)
3. Sampling coordinates
4. Environmental contrasts among habitats
5. Taxonomic summaries
6. Family–environment associations
7. Taxonomic resolution and community overview
8. Taxonomic composition
9. Alpha diversity and Hill number profiles
10. Ordination (PCoA, NMDS, RDA) and environmental models
11. Genus association networks (SpiecEasi) and visualisation
12. Network robustness simulations
13. Within-module and among-module node roles
14. Mantel tests and environmental correlations
15. Network topology summary
16. Supplementary taxonomic and guild tables

## How to run

1. Install R and the packages listed below. Running `source("R/packages.R")` once installs anything that is missing.
2. Place the input file at `data/fungi.RData` (see `data/README.md` for the expected format).
3. Open `fungi-analysis.Rproj` in RStudio, or otherwise start R with the repository root as the working directory.
4. Run `source("R/Fungi_Paper_Script.R")`.

All paths are relative to the repository root. Figures (PDF and PNG) and tables (CSV and XLSX) are written to `outputs/`, which the script creates if needed and which is not tracked by git. The analysis settings, including the random seed (`123`), number of permutations (`999`) and the switch for the optional ANOVA comparison, are set in the configuration block at the top of the script. At the end of a run, the script writes the R session details, including package versions, to `outputs/session_info.txt`.

## R packages

| Source | Packages |
|---|---|
| CRAN | dplyr, tidyr, purrr, tibble, ggplot2, rstatix, ggpubr, patchwork, Polychrome, Hmisc, ggvenn, circlize, writexl, vegan, ggforce, igraph, ggraph, tidygraph, ggnewscale, scales, ggrepel, cowplot, FSA, rcompanion |
| CRAN (only when `run_anova_comparison <- TRUE`) | multcompView, car |
| Bioconductor | phyloseq, ComplexHeatmap |
| GitHub | SpiecEasi, linkET |

The base packages `grid` and `stats` are also used. The R and package versions used for the manuscript are recorded in `session_info.txt` [ADD THIS FILE, COPIED FROM `outputs/` AFTER THE FINAL RUN].

## Citation

If you use this code, please cite the manuscript: [CITATION TO BE ADDED ON PUBLICATION].

## License

[LICENSE TO BE ADDED]

## Contact

Rapheal Wangalwa, [Department of Biology, Faculty of Science, Mbarara University of Science and Technology, Uganda; Email: wangarapha@must.ac.ug]

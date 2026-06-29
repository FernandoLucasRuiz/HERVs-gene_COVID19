pkgs <- c(
    "tidyverse", "knitr", "readxl", "DESeq2", "gtsummary", "tximport",
    "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "circlize",
    "tidyheatmaps", "viridisLite", "clusterProfiler", "ReactomePA",
    "DOSE", "biomaRt", "GWENA", "ggsankey", "ggpubr", "ggrepel",
    "igraph", "scales"
)

invisible(lapply(pkgs, library, character.only = TRUE))

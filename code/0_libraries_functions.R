
library(tidyverse)
library(knitr)
library(readxl)
library(DESeq2) 

library(gtsummary)

library(tximport) 

library(TxDb.Hsapiens.UCSC.hg38.knownGene)

library(org.Hs.eg.db)

library(circlize)

library(tidyheatmaps)

library(viridisLite)

library(clusterProfiler)
library(ReactomePA)
library(DOSE)

library(biomaRt)

library(GWENA)

library(ggsankey)
library(ggpubr)
library(ggrepel)

library(igraph)
library(scales)

library(Seurat)
library(AnnotationDbi)
library(SingleCellExperiment)
library(zellkonverter)

library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(ComplexUpset)
library(corrplot)

colores <- c(
    "#2A9D8F", "#264653",  "#8AB17D", "#E9C46A", "#F4A261",
    "#E76F51", "#F28482", "#9B5DE5", "#4E8098", "#577590",
    "#43AA8B", "#90BE6D", "#F8961E", "#F3722C", "#277DA1",
    "#6D597A", "#B56576", "#EAAC8B", "#84A59D", "#CDB4DB",
    "#A8DADC"
    )

color_control = colores[9]
color_covid= colores[17]

plot_enriquecimientos <- function(df, GeneRatio_num, Description, pvalue, titulo = NULL) {
    
    p <- ggplot(df, 
                aes(
                    x = {{ GeneRatio_num }}, 
                    y = {{ Description }}, 
                    color = -log10({{ pvalue }} ))) + 
    geom_segment(aes(x=0, xend={{ GeneRatio_num }}, y={{ Description }}, yend={{ Description }}), 
                 color="skyblue",
                 size = 1) +
    geom_point(size = 4) +
    scale_size(range = c(4, 7)) +
    scale_color_gradient2(low = "#2A9D8F", mid = "white", high = "#A52A2A") +
    facet_grid(
        rows = vars(db),
        labeller = labeller(db = lab_names),
        scales = "free_y",
        space = "free_y"
    ) +
    labs(
        y = NULL,
        x = "Gene ratio",
        color = expression(-log[10]("pvalue")),
        title = titulo
    ) +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        strip.text.y = element_text(
            size = 10,
            face = "bold",
            angle = 0,
            hjust = 0
            ),
        strip.background = element_blank(),
        plot.margin = margin(20,20,20,20)
    )
    return(p)
} 

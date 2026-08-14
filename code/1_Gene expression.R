
ddsSalmon_filtered_GC <- readRDS("../data/dds_DEG_filtered_GC.rds")
ddsSalmon_filtered_GC <- DESeq(ddsSalmon_filtered_GC)
res_DEG_uni_GC <- results(ddsSalmon_filtered_GC, contrast = c("condition", "COVID", "Control"))

```

## Annotation of results table

```{r}
Gene_ann_GC <- AnnotationDbi::select(org.Hs.eg.db,
                 keys = rownames(res_DEG_uni_GC),
                 keytype = "ENTREZID",
                 columns="SYMBOL")

res_DEG_uni_ann_GC <- merge(Gene_ann_GC,as.data.frame(res_DEG_uni_GC),
                      by.x=1,
                      by.y=0,
                      all.x=FALSE,
                      all.y=TRUE)

res_DEG_uni_ann_GC <- res_DEG_uni_ann_GC |>
    mutate(
        expression= case_when(
            log2FoldChange <= -1 & padj <= 0.05 ~ "significative downregulated",
            log2FoldChange <= -1 & padj > 0.05 ~ "downregulated",
            log2FoldChange >= 1 & padj < 0.05 ~ "significative upregulated",
            log2FoldChange >= 1 & padj > 0.05 ~ "downregulated",
            log2FoldChange < 1 & log2FoldChange > -1 & padj < 0.05 ~ "significative",
            log2FoldChange < 1 & log2FoldChange > -1 & padj > 0.05 ~ "nothing"
        )
    )


genes_DE_anotados <- res_DEG_uni_ann_GC

rm(res_DEG_uni_ann_GC)

```

## Tabla1

```{r}
genes_DE_anotados |>
    filter(abs(log2FoldChange) >=1 & padj <= 0.05) |>
    dplyr::select(SYMBOL, log2FoldChange, padj) |> write.csv("../../../../Artículos/HERVs/version0_no20E_RINadj/tablas/SupplementaryTable1.csv")
```

## Housekeeping

```{r}
hk_candidates <- intersect(
  genes_DE_anotados$SYMBOL,
  c(
      "ACTB", "APRT", "EEF1A1", "UBC", "HPRT1", "SDHA", "TUBB", "PPIA","TBP","PGK1","GUSB","PSMB4","VCP","POLR2A", "B2M", "RPLP1", "RPLP2", "TFRC", "ALAS1", "PUM1",
        "YWHAZ", "IPO8",
  "UBE2D2", "C1orf43", "EMC7"
  )
)
```

```{r, fig.width=4, fig.height=4}
Genes_mcount_vst <- vst(assay(ddsSalmon_filtered_GC))
genes_entrezid_hk <- genes_DE_anotados |>
    filter(SYMBOL %in% hk_candidates) |> pull(ENTREZID)

mat <- Genes_mcount_vst[rownames(Genes_mcount_vst) %in% genes_entrezid_hk,]

df <- mat |>
    t() |> 
    as.data.frame() |> 
    rownames_to_column("sample") |> 
    mutate(across(where(is.factor), as.character)) |>  # Convierte factores a caracteres
    pivot_longer(cols = -sample) |>
    dplyr::left_join(patients_data_all_cov_RIN, by = c("sample" = "expediente")) |>
    dplyr::left_join(genes_DE_anotados, by = c("name" = "ENTREZID"))


ggboxplot(df, x = "condition", y = "value",
          fill = "condition") +
    geom_jitter(alpha=0.5, width = 0.15) +
    stat_compare_means(label = "p.format", hjust = 0.5,label.x = 1.5) +
    scale_fill_manual(values = c(color_control, color_covid)) +
    theme_minimal() +
    labs(
        x = NULL,
        y= "Normalized HK genes expression"
    ) +
    theme(
        axis.text.x = element_text(face="bold"),
        plot.margin = margin(20,20,20,20), 
        legend.position = "none"
    )

```

```{r, fig.width=10, fig.height=8}
ggboxplot(df, x = "condition", y = "value", fill = "condition", add = "jitter") +
    scale_fill_manual(values = c(color_control, color_covid)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    facet_wrap(~SYMBOL, scales = "free") +
    theme_bw() +
    stat_compare_means(label = "p.signif", label.x = 1.5,hjust = 0.5)  +
    labs(
        x= NULL,
        y= "Normalized gene expression"
    ) +
    theme(
        legend.position = "none"
    )
```

## heatmap

```{r, fig.width=6, fig.height=4}
vsd <- vst(ddsSalmon_filtered_GC, blind = FALSE)
mat <- assay(vsd)

df <- as.data.frame(mat) |>
    rownames_to_column("genes") |>
    pivot_longer(-genes, names_to = "sample", values_to = "expression")

metadata <- as.data.frame(colData(ddsSalmon_filtered_GC)) |>
    rownames_to_column("sample")

df <- df |>
    left_join(metadata, by = "sample") 

# top50_COVID_vs_control <- HERVs_DE_anotados |>
#     arrange(desc(log2FoldChange)) |>
#     slice_head(n = 50) |>
#     pull(HERVs)

tidyheatmap(
    df,
    rows = genes,
    columns = sample,
    values = expression,
    scale = "row", 
    color_legend_min = -1,
    color_legend_max = 1,
    colors = c("#2C7BB6","#F7F7F7","#FDAE61"),
    color_legend_n = 100,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = condition,
    clustering_method = "ward.D2",
    cutree_cols = 2,
    show_rownames = F,
    annotation_colors = list(
        condition = c(Control= color_control, COVID =color_covid)
        )
    ) 
```

## umap

```{r, fig.width=5, fig.height=3}
vsd <- vst(ddsSalmon_filtered_GC, blind = FALSE)
mat <- t(assay(vsd))

set.seed(1234)
um <- uwot::umap(mat)

metadata <- as.data.frame(colData(ddsSalmon_filtered_GC)) |>
  rownames_to_column("sample")

um_df <- um |> as.data.frame() |>
    rownames_to_column("sample") |>
    left_join(metadata, by = "sample") 

centroids <- um_df %>%
    group_by(condition) %>%
    summarize(V1 = median(V1, na.rm = TRUE),
              V2 = median(V2, na.rm = TRUE))

labels_manual <- data.frame(
    conditio = c("Control", "COVID"),
    X1 = c(-7, 5.5),   # posiciones X elegidas por ti
    X2 = c(-3, -4)   # posiciones Y elegidas por ti
    )

um_df|>
    ggplot(aes(x=V1, y=V2, color=condition)) +
    geom_point(size=4)+
    labs(
        x="UMAP1",
        y="UMAP2",
        color=NULL
    ) +
    stat_ellipse(data = um_df, aes(x = V1, y = V2, fill = condition),
                 level = 0.95,
                 show.legend = F,
                 geom = "polygon", color = NA,
                 alpha = 0.2) +
    scale_fill_manual(values = c(colores[9], colores[17])) +
    scale_color_manual(values = c(colores[9], colores[17])) +
    theme_minimal() +
    theme(
        legend.position = "bottom", 
        plot.margin = margin(20, 20, 20, 20)
    )
```


## Volcano

```{r, fig.width=6, fig.height=5}
most_fc <- genes_DE_anotados |>
    filter(abs(log2FoldChange) >= 1.5 & padj <= 0.05) |>
    arrange(desc(abs(log2FoldChange))) |>
    slice_head(n = 10) |>
    pull(SYMBOL)
most_padj <- genes_DE_anotados |>
    filter(abs(log2FoldChange) >= 1.5 & padj <= 0.05) |>
    arrange(padj) |>
    slice_head(n = 10) |>
    pull(SYMBOL)

labels <- c(most_fc, most_padj)
tab <- table(genes_DE_anotados$expression)

genes_DE_anotados |>
    ggplot(aes(x=log2FoldChange, y=-log10(padj), color = expression)) +
    geom_point(
        data = \(x) dplyr::filter(x, expression %in% c("significative upregulated", "significative downregulated")),
        aes(color = expression),
        alpha = 1,
        size = 3
    ) +
      geom_point(
        data = \(x) dplyr::filter(x, !expression %in% c("significative upregulated", "significative downregulated")),
        aes(color = expression),
        alpha = 0.5,
        size = 1
      ) +
    # ggrepel::geom_text_repel(
    #     data = ~ dplyr::filter(.x, SYMBOL %in% labels),
    #     aes(label = SYMBOL),
    #     size = 3,
    #     color = "black",
    #   
    #     force = 5,              # más repulsión entre etiquetas
    #     force_pull = 1,         # mantiene cierta conexión con el punto
    #     box.padding = 0.6,      # espacio alrededor del texto
    #     point.padding = 1.2,    # separación respecto al punto
    #     max.overlaps = Inf,     # no descarta etiquetas
    #   
    #     min.segment.length = 0, # siempre dibuja línea
    #     segment.color = "grey50",
    #     segment.size = 0.3,
    #   
    #     max.time = 2,           # más tiempo para optimizar posiciones
    #     max.iter = 20000        # más iteraciones → mejor colocación
    #     ) +
    annotate(
        "text",
        x = 3.5,
        y = -3,
        label = paste0(tab[5], " Upregulated genes "),
        size = 4,
        vjust = 0.5,
        hjust = 0.5, 
        fontface = "bold",
        color=color_covid
    ) +
    annotate(
        "text",
        x = -3.2,
        y = -3,
        label = paste0(tab[4], " Downregulated genes "),
        size = 4,
        vjust = 0.5,
        hjust = 0.5, 
        fontface = "bold",
        color=color_control
    ) +
    scale_color_manual(values = c(colores[21], colores[21], colores[21], color_control, color_covid)) +
    theme_minimal() +
    labs(
        y=expression(-log[10]("adj. p")),
        x=expression(log[2]("FC")),
        color = NULL
    ) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dotted") +
    geom_vline(xintercept = c(-1, 1), color = "red", linetype = "dotted") +
    xlim(c(-7,7)) +
    ylim(c(-4,30)) +

    theme(
        legend.position = "bottom", 
        margins = margin(20, 20, 20, 20)
    )
```

### Housekeeping

```{r, fig.width=5, fig.height=5}
labels <- hk_candidates

genes_DE_anotados |>
  dplyr::mutate(
    housekeepings = SYMBOL %in% hk_candidates
  ) |>
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  
  geom_point(
    data = \(x) dplyr::filter(x, !housekeepings),
    color = "gray",
    alpha = 0.1,
    size = 1
  ) +
  
  geom_point(
    data = \(x) dplyr::filter(x, housekeepings),
    color = "red",
    alpha = 1,
    size = 2.5
  ) +
  
  ggrepel::geom_text_repel(
    data = \(x) dplyr::filter(x, SYMBOL %in% labels),
    aes(label = SYMBOL),
    size = 3,
    color = "black",
  
    force = 5,              # más repulsión entre etiquetas
    force_pull = 1,         # mantiene cierta conexión con el punto
    box.padding = 0.6,      # espacio alrededor del texto
    point.padding = 0.5,    # separación respecto al punto
    max.overlaps = Inf,     # no descarta etiquetas
  
    min.segment.length = 0, # siempre dibuja línea
    segment.color = "grey50",
    segment.size = 0.3,
  
    max.time = 2,           # más tiempo para optimizar posiciones
    max.iter = 20000        # más iteraciones → mejor colocación
    ) +
  
  theme_minimal() +
  labs(
    y = expression(-log[10]("adj. p")),
    x = expression(log[2]("FC"))
  ) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dotted") +
  geom_vline(xintercept = c(-1, 1), color = "red", linetype = "dotted") +
      xlim(c(-7,7)) +
    ylim(c(-4,30)) +
  theme(
    legend.position = "none",
    plot.margin = margin(20, 20, 20, 20)
  )
```

## Functional enrichment

### Normal

```{r, eval=FALSE}
enriquecimientos <- list()

lista_enriquecer.df <- genes_DE_anotados |>
    filter(abs(log2FoldChange) >= 1.5 & padj <= 0.05)
    
ego <- enrichGO(gene = unique(lista_enriquecer.df$ENTREZID),
                        keyType       = "ENTREZID",
                        OrgDb         = org.Hs.eg.db,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.1,
                        readable      = TRUE)
    
enriquecimientos[["GO"]] <- ego
    
kk <- enrichKEGG(gene         = unique(lista_enriquecer.df$ENTREZID), 
                     organism     = 'hsa',
                     pvalueCutoff = 0.1)

enriquecimientos[["KEGG"]] <- kk

wp <- enrichWP(gene = unique(lista_enriquecer.df$ENTREZID),
         organism = "Homo sapiens") 
    
enriquecimientos[["WikiPathways"]] <- wp

x <- enrichPathway(gene=unique(lista_enriquecer.df$ENTREZID), 
                   pvalueCutoff = 0.1, 
                   readable=TRUE)

enriquecimientos[["Reactome"]] <- x

# x <- enrichDO(gene          = unique(lista_enriquecer.df$ENTREZID),
#               ont           = "HDO",
#               pvalueCutoff  = 0.1,
#               pAdjustMethod = "BH",
#               universe      = names(lista_enriquecer.df$ENTREZID),
#               minGSSize     = 10,
#               maxGSSize     = 500,
#               qvalueCutoff  = 0.05,
#               readable      = FALSE)
# 
# enriquecimientos[["DOSE_Disease"]] <- x
# 
# x <- enrichDO(gene          = unique(lista_enriquecer.df$ENTREZID),
#               ont           = "HPO",
#               pvalueCutoff  = 0.1,
#               pAdjustMethod = "BH",
#               universe      = names(lista_enriquecer.df$ENTREZID),
#               minGSSize     = 10,
#               maxGSSize     = 500,
#               qvalueCutoff  = 0.05,
#               readable      = FALSE)
# 
# enriquecimientos[["DOSE_Phenotype"]] <- x

saveRDS(enriquecimientos, "data/no_20E_RINadj/enriquecimientos_DEGs.RDS")
```

```{r}
enriquecimientos <- readRDS("data/no_20E_RINadj/enriquecimientos_DEGs.RDS")
```




```{r}
enriquecimientos_filtered_df <- data.frame()
enriquecimientos_filtrados <- enriquecimientos

for (i in names(enriquecimientos)) {
    
    df_filtered <- enriquecimientos[[i]]@result #|>
        # filter(grepl('*chrom*', Description) | grepl('*methyl*', Description) | grepl('*stimulus*', Description)) 
    
    enriquecimientos_filtered_df <- df_filtered |> 
        dplyr::select(Description, Count, GeneRatio, p.adjust, pvalue, geneID) |>
        mutate(db = i) |> 
        rbind(enriquecimientos_filtered_df)
    
    
}

enriquecimientos_filtered_df <- enriquecimientos_filtered_df |> filter(p.adjust <= 0.05)
```

```{r}
terminos_cromatina <- c(
  "Chromatin modifications during the maternal to zygotic transition (MZT)",
  "DNA methylation",
  "chromosome condensation",
  "structural constituent of chromatin"
)

terminos_antiviral <- c(
  "Interferon alpha/beta signaling",
  "Type I interferon induction and signaling during SARS CoV 2 infection",
  "SARS CoV 2 innate immunity evasion and cell specific immune response",
  "response to type I interferon",
  "cellular response to type I interferon",
  "type I interferon-mediated signaling pathway",
  "Interferon Signaling",
  "defense response to virus",
  "response to virus",
  "cellular response to dsRNA",
  "interferon-mediated signaling pathway",
  "positive regulation of type I interferon production",
  "positive regulation of interferon-beta production",
  "viral process",
  "viral genome replication",
  "regulation of viral process",
  "negative regulation of viral process",
  "negative regulation of viral genome replication"
)

terminos_inflamacion <- c(
  "positive regulation of cytokine production",
  "positive regulation of cytokine production involved in immune response",
  "response to lipopolysaccharide",
  "response to molecule of bacterial origin",
  "response to chemokine",
  "cellular response to chemokine",
  "chemokine-mediated signaling pathway",
  "leukocyte migration",
  "leukocyte chemotaxis",
  "leukocyte migration involved in inflammatory response",
  "myeloid leukocyte migration",
  "mononuclear cell migration",
  "macrophage migration",
  "neutrophil migration",
  "neutrophil chemotaxis",
  "complement activation",
  "complement activation, classical pathway",
  "positive regulation of tumor necrosis factor production",
  "positive regulation of tumor necrosis factor superfamily cytokine production",
  "interleukin-1 beta production",
  "regulation of interleukin-1 beta production"
)

terminos_ecm_fibrosis <- c(
  "Extracellular matrix organization",
  "external encapsulating structure organization",
  "extracellular matrix assembly",
  "extracellular matrix constituent secretion",
  "ECM proteoglycans",
  "ECM-receptor interaction",
  "Collagen formation",
  "Collagen biosynthesis and modifying enzymes",
  "collagen biosynthetic process",
  "collagen metabolic process",
  "collagen fibril organization",
  "Collagen chain trimerization",
  "Collagen degradation",
  "collagen catabolic process",
  "collagen-containing extracellular matrix",
  "Lung fibrosis",
  "tissue remodeling",
  "wound healing",
  "fibroblast proliferation",
  "epithelial to mesenchymal transition"
)
```

```{r,fig.width=8, fig.height=8}
lab_names <- c(
  GO = "Gene Ontology",
  KEGG = "KEGG",
  Reactome = "Reactome",
  # DOSE_Disease = "Disease Ontology",
  # DOSE_Phenotype = "Phenotype Ontology",
  WikiPathways = "WikiPathways"
)

df_enriquecido <- enriquecimientos_filtered_df |>
  group_by(db) |>
    mutate(
    db = factor(
      db,
      levels = c(
        "GO",
        "KEGG",
        "Reactome",
        # "DOSE_Disease",
        # "DOSE_Phenotype",
        "WikiPathways"
      )
    )
  )|>
    mutate(log10_p = -log10(pvalue)) |>
    mutate(
        Description = if_else(
            str_length(Description) > 50,
            str_c(str_sub(Description, 1, 50), "…"),
            Description
            )
    )|>
    mutate(
        GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) /
            as.numeric(sub(".*/", "", GeneRatio))
    )

```

```{r,fig.width=8, fig.height=12}
# epigenetica_cromatina, regulacion_transcripcional, senalizacion_pro_transcripcional, arquitectura_cromatina, inflamasoma_core, inflamacion_directa, inflamacion_enfermedad


plot_df <- df_enriquecido |>
    filter(Description %in% c(epigenetica_cromatina,
                              regulacion_transcripcional,
                              senalizacion_pro_transcripcional,
                              arquitectura_cromatina)) |>
     mutate(
        Description = fct_reorder(Description, pvalue, .desc = TRUE)
    )

plot_enriquecimientos(df = plot_df,
                      GeneRatio_num = GeneRatio_num ,
                      Description = Description, 
                      pvalue = pvalue, 
                      titulo = "Chromatin & trancriptional terms") +
    labs(
        size = NULL
    ) +
    theme(
        axis.text.y = element_text(size=12), 
        title = element_text(size = 16, face = "bold"),
        strip.text.y = element_text(
            size = 14,   # aumenta el tamaño
            face = "bold"
        )
    ) 

```

```{r,fig.width=8, fig.height=12}
# epigenetica_cromatina, regulacion_transcripcional, senalizacion_pro_transcripcional, arquitectura_cromatina, inflamasoma_core, inflamacion_directa, inflamacion_enfermedad

plot_df <- df_enriquecido |>
    filter(Description %in% c(inflamasoma_core, 
                              inflamacion_directa)) |>
     mutate(
        Description = fct_reorder(Description, pvalue, .desc = TRUE)
    )

plot_enriquecimientos(df = plot_df,
                      GeneRatio_num = GeneRatio_num ,
                      Description = Description, 
                      pvalue = pvalue, 
                      titulo = "Inflammation") +
    theme(
        axis.text.y = element_text(size=12), 
        title = element_text(size = 16, face = "bold"),
        strip.text.y = element_text(
            size = 14,   # aumenta el tamaño
            face = "bold"
        )
    ) 

```

```{r, fig.width=10, fig.height=4}
genes_unicos_cromatina <- enriquecimientos_filtered_df %>%
    filter(Description %in% c(epigenetica_cromatina,
                              regulacion_transcripcional,
                              senalizacion_pro_transcripcional,
                              arquitectura_cromatina)) |>
    pull(geneID) %>%
    strsplit("/") %>%
    unlist() %>%
    unique()
genes_unicos_cromatina <- genes_DE_anotados |>
    filter(SYMBOL %in% genes_unicos_cromatina | ENTREZID %in% genes_unicos_cromatina) 

genes_unicos_inflamacion <- df_enriquecido %>%
    filter(Description %in% c(inflamasoma_core, 
                              inflamacion_directa)) |>
    pull(geneID) %>%
    strsplit("/") %>%
    unlist() %>%
    unique()
genes_unicos_inflamacion <- genes_DE_anotados |>
    filter(SYMBOL %in% genes_unicos_inflamacion | ENTREZID %in% genes_unicos_inflamacion)

vsd <- vst(ddsSalmon_filtered_GC, blind = FALSE)
mat <- assay(vsd) |> as.data.frame()

mat |>
  rownames_to_column("genes") |>
  filter(genes %in% c(genes_unicos_cromatina$ENTREZID, genes_unicos_inflamacion$ENTREZID)) |>
  pivot_longer(-genes) |>
    mutate(termino = ifelse(genes %in% genes_unicos_cromatina$ENTREZID, 
                            "Chromatin terms", 
                            "Inflammation terms")) |>
  left_join(patients_data_all_cov, by = c("name" = "expediente")) |>
    ggplot(aes(x=condition, y = value, fill = condition)) +
    geom_boxplot() +
    geom_point(
        aes(x=condition, y = value, fill = condition),
        shape =21,
        alpha=0.5,
        position = position_jitter(width = 0.2)) +
    scale_fill_manual(values = c(color_control, color_covid)) +
    scale_color_manual(values = c(color_control, color_covid)) +
    stat_compare_means(
        label = "p.signif", 
        comparisons = list(c("Control", "COVID")),
        label.x = 1.5, 
        hjust = 0.5
        ) +
    facet_wrap(~termino, scales = "free_y") +
    theme_minimal() +
    labs(
        x=NULL, 
        y="Normalized gene expression"
    ) +
    theme(
        strip.text  = element_text(size=16, face="bold"),
        axis.text.x = element_text(size = 10, face="bold"),
        legend.position = "none",
        plot.margin = margin(20,20,20,20)
    ) 

```

#### Tabla enriquecimientos genes

```{r}
# df_enriquecido |>
#     ungroup() |>
#     filter(Description %in% c(inflamasoma_core, 
#                               inflamacion_directa, 
#                               epigenetica_cromatina,
#                               regulacion_transcripcional,
#                               senalizacion_pro_transcripcional,
#                               arquitectura_cromatina)) |>
#     dplyr::select(Description, GeneRatio, p.adjust, db) |> 
#     write.csv("../../../Artículos/HERVs/version0/tablas/SupplementaryTable2.csv")
```

### gsea

```{r}
library(dplyr)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

#------------------------------------------------------------
# 1. Crear la lista ordenada de genes
#------------------------------------------------------------

lista_GSEA.df <- genes_DE_anotados |>
    filter(
        !is.na(ENTREZID),
        !is.na(stat),
        is.finite(stat)
    ) |>
    mutate(
        ENTREZID = as.character(ENTREZID)
    ) |>
    group_by(ENTREZID) |>
  # Si varios genes/Ensembl corresponden al mismo ENTREZID,
  # conservar el que tenga el estadístico más extremo
    slice_max(
        order_by = abs(stat),
        n = 1,
        with_ties = FALSE
    ) |>
    ungroup()

geneList <- lista_GSEA.df$stat
names(geneList) <- lista_GSEA.df$ENTREZID

# GSEA necesita la lista ordenada de mayor a menor
geneList <- sort(geneList, decreasing = TRUE)

set.seed(123)

enriquecimientos_GSEA <- list()

#------------------------------------------------------------
# 2. Gene Ontology
#------------------------------------------------------------

gsea_GO <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  exponent      = 1,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.1,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

enriquecimientos_GSEA[["GO"]] <- gsea_GO


#------------------------------------------------------------
# 3. KEGG
#------------------------------------------------------------

gsea_KEGG <- gseKEGG(
  geneList      = geneList,
  organism      = "hsa",
  keyType       = "ncbi-geneid",
  exponent      = 1,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.1,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

enriquecimientos_GSEA[["KEGG"]] <- gsea_KEGG


#------------------------------------------------------------
# 4. WikiPathways
#------------------------------------------------------------

gsea_WP <- gseWP(
  geneList      = geneList,
  organism      = "Homo sapiens",
  exponent      = 1,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.1,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

enriquecimientos_GSEA[["WikiPathways"]] <- gsea_WP


#------------------------------------------------------------
# 5. Reactome
#------------------------------------------------------------

gsea_Reactome <- gsePathway(
  geneList      = geneList,
  organism      = "human",
  exponent      = 1,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.1,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

enriquecimientos_GSEA[["Reactome"]] <- gsea_Reactome

saveRDS(enriquecimientos_GSEA, "../data/no_20E_RINadj/enriquecimientos_GSEA_DEGs.RDS")
```

```{r}
enriquecimientos <- readRDS("../data/no_20E_RINadj/enriquecimientos_GSEA_DEGs.RDS")
```

```{r}
enriquecimientos_filtered_df <- data.frame()
enriquecimientos_filtrados <- enriquecimientos

for (i in names(enriquecimientos)) {
    
    df_filtered <- enriquecimientos[[i]]@result #|>
        # filter(grepl('*chrom*', Description) | grepl('*methyl*', Description) | grepl('*stimulus*', Description)) 
    
    enriquecimientos_filtered_df <- df_filtered |> 
        dplyr::select(Description, NES, p.adjust, pvalue) |>
        mutate(db = i) |> 
        rbind(enriquecimientos_filtered_df)
    
    
}

enriquecimientos_filtered_df <- enriquecimientos_filtered_df |> filter(p.adjust <= 0.05)
```

```{r}
gsea_cromatina_epigenetica_core <- c(
  "DNA methylation",
  "HDACs deacetylate histones",
  "PRC2 methylates histones and DNA",
  "RMTs methylate histone arginines",
  "Chromatin modifications during the maternal to zygotic transition (MZT)",
  "structural constituent of chromatin",
  "nucleosome",
  "Nucleosome assembly",
  "Deposition of new CENPA-containing nucleosomes at the centromere",
  "Orc1 removal from chromatin",
  "Chromosome Maintenance",
  "Condensation of Prophase Chromosomes",
  "Condensation of Prometaphase Chromosomes",
  "chromosome segregation",
  "nuclear chromosome segregation",
  "chromosome separation",
  "regulation of chromosome separation",
  "negative regulation of chromosome organization",
  "protein localization to condensed chromosome",
  "Telomere Maintenance",
  "Packaging Of Telomere Ends",
  "Inhibition of DNA recombination at telomere"
)

gsea_antiviral_interferon_core <- c(
  "Interferon alpha/beta signaling",
  "Interferon Signaling",
  "Type II interferon signaling",
  "DDX58/IFIH1-mediated induction of interferon-alpha/beta",
  "Antiviral mechanism by IFN-stimulated genes",
  "ISG15 antiviral mechanism",
  "response to interferon-beta",
  "response to type I interferon",
  "cellular response to type I interferon",
  "type I interferon-mediated signaling pathway",
  "interferon-mediated signaling pathway",
  "positive regulation of type I interferon production",
  "positive regulation of interferon-beta production",
  "defense response to virus",
  "response to virus",
  "antiviral innate immune response",
  "regulation of viral process",
  "negative regulation of viral process",
  "viral process",
  "viral genome replication",
  "regulation of viral genome replication",
  "negative regulation of viral genome replication",
  "regulation of viral life cycle",
  "negative regulation of viral life cycle",
  "viral life cycle"
)

gsea_viral_retroviral_herv <- c(
  "Host Interactions of HIV factors",
  "HIV Infection",
  "Vpu mediated degradation of CD4",
  "Vif-mediated degradation of APOBEC3G",
  "Interactions of Rev with host cellular proteins",
  "Rev-mediated nuclear export of HIV RNA",
  "Early Phase of HIV Life Cycle",
  "HIV Life Cycle",
  "The role of Nef in HIV-1 replication and disease pathogenesis",
  "Viral mRNA Translation",
  "Viral Messenger RNA Synthesis",
  "Export of Viral Ribonucleoproteins from Nucleus",
  "Influenza Infection",
  "Influenza Viral RNA Transcription and Replication",
  "SARS-CoV-2 Infection",
  "SARS-CoV-2-host interactions",
  "SARS-CoV-2 modulates host translation machinery",
  "SARS-CoV-2 activates/modulates innate and adaptive immune responses",
  "Coronavirus disease",
  "Network map of SARS CoV 2 signaling"
)

gsea_inflamacion_persistente_core <- c(
  "Dectin-1 mediated noncanonical NF-kB signaling",
  "TNFR2 non-canonical NF-kB pathway",
  "Activation of NF-kappaB in B cells",
  "CLEC7A (Dectin-1) signaling",
  "Interleukin-1 signaling",
  "Interleukin-1 family signaling",
  "Signaling by Interleukins",
  "Defective pyroptosis",
  "Initial triggering of complement",
  "Downstream signaling events of B Cell Receptor (BCR)",
  "Downstream TCR signaling"
)

gsea_estres_oxidativo_hipoxia <- c(
  "Nuclear events mediated by NFE2L2",
  "KEAP1-NFE2L2 pathway",
  "GSK3B and BTRC:CUL1-mediated-degradation of NFE2L2",
  "NRF2 pathway",
  "Transcriptional activation by NRF2 in response to phytochemicals",
  "Oxygen-dependent proline hydroxylation of Hypoxia-inducible Factor Alpha",
  "Cellular response to hypoxia",
  "Chemical carcinogenesis - reactive oxygen species",
  "Oxidative Stress Induced Senescence"
)

gsea_ecm_fibrosis_remodelado <- c(
  "Extracellular matrix organization",
  "extracellular matrix organization",
  "extracellular matrix",
  "collagen fibril organization",
  "collagen-containing extracellular matrix",
  "extracellular matrix structural constituent",
  "extracellular matrix structural constituent conferring tensile strength",
  "extracellular matrix structural constituent conferring compression resistance",
  "Collagen formation",
  "Collagen biosynthesis and modifying enzymes",
  "Collagen degradation",
  "Degradation of the extracellular matrix",
  "Assembly of collagen fibrils and other multimeric structures",
  "Collagen chain trimerization",
  "collagen metabolic process",
  "collagen catabolic process",
  "collagen binding",
  "collagen trimer",
  "fibronectin binding",
  "fibroblast proliferation",
  "ECM proteoglycans",
  "Integrin cell surface interactions",
  "Elastic fibre formation",
  "Activation of Matrix Metalloproteinases",
  "Matrix metalloproteinases",
  "Type I collagen synthesis in the context of osteogenesis imperfecta",
  "miRNA targets in ECM and membrane receptors",
  "miR 509 3p alteration of YAP1 ECM axis",
  "Metabolic pathways of fibroblasts",
  "Platelet Adhesion to exposed collagen"
)
```

```{r}
gsea_df <- enriquecimientos_filtered_df

term_categories_gsea <- tibble::tibble(
  Description = c(
    gsea_cromatina_epigenetica_core,
    gsea_antiviral_interferon_core,
    gsea_viral_retroviral_herv,
    gsea_inflamacion_persistente_core,
    gsea_estres_oxidativo_hipoxia,
    gsea_ecm_fibrosis_remodelado
  ),
  Category = c(
    rep("Chromatin / epigenetic regulation", length(gsea_cromatina_epigenetica_core)),
    rep("Antiviral / interferon response", length(gsea_antiviral_interferon_core)),
    rep("Viral / retroviral-like programs", length(gsea_viral_retroviral_herv)),
    rep("Persistent inflammation", length(gsea_inflamacion_persistente_core)),
    rep("Oxidative stress", length(gsea_estres_oxidativo_hipoxia)),
    rep("ECM remodeling / fibrosis", length(gsea_ecm_fibrosis_remodelado))
  )
) |>
  mutate(Description_clean = str_squish(as.character(Description))) |>
  distinct(Description_clean, Category)

plot_df_gsea <- gsea_df |>
  mutate(
    Description_clean = str_squish(as.character(Description)),
    NES = as.numeric(NES),
    p.adjust = as.numeric(p.adjust),
    pvalue = as.numeric(pvalue),
    minus_log10_padj = -log10(pmax(p.adjust, 1e-300)),
    Direction = ifelse(NES > 0, "Upregulated genes", "Downregulated genes")
  ) |>
  inner_join(term_categories_gsea, by = "Description_clean") |>
    mutate(
        Category = factor(Category, levels = c("Chromatin / epigenetic regulation", "Persistent inflammation", "Antiviral / interferon response","Viral / retroviral-like programs","Oxidative stress","ECM remodeling / fibrosis"
))
    )
  

```

```{r,fig.width=10, fig.height=12}
plot_df_gsea |>
    group_by(Category) |>
  slice_min(order_by = p.adjust, n = 8, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    Description_plot = fct_reorder(Description_clean, NES)
  ) |>
    ggplot(aes(x = NES, y = Description_plot)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_segment(
    aes(x = 0, xend = NES, yend = Description_plot),
    linewidth = 0.6,
    colour = "grey75"
  ) +
  geom_point(
    aes(size = minus_log10_padj, fill = minus_log10_padj),
    shape = 21,
    colour = "black"
  ) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(
    low = "grey85",
    high = "#A51E22"
  ) +
  labs(
    x = "Normalized enrichment score (NES)",
    y = NULL,
    fill = expression(-log[10](FDR)),
    size = expression(-log[10](FDR))
  ) +
  theme_bw(base_size = 12) +
  theme(
       axis.text.y = element_text(size = 10),
    strip.text.y = element_text(size= 12,angle = 0, face = "bold"),
    strip.background = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5)
   
  )
```

```{r, fig.width=10, fig.height=3}
plot_df_gsea |>
    group_by(Category) |>
    summarize(Count = n()) |>
    mutate(Category = fct_reorder(Category, Count)) |>
    ggplot(aes(y = Category, x = Count)) +
    geom_col(fill = colores[10]) +
    theme_minimal() +
    labs(
        y=NULL,
        x= "Number of enriched terms"
    )+
    theme(
        axis.text.y = element_text(size= 12, face="bold"),
        axis.title.x  = element_text(size= 12, face="bold"),
    )
```


```{r,fig.width=5, fig.height=3}
library(enrichplot)

# Buscar el ID real del término dentro del objeto GSEA
id <- enriquecimientos$GO@result |>
  filter(Description == "structural constituent of chromatin") |>
  pull(ID)

gseaplot(
    enriquecimientos$GO,
    geneSetID = id,
    color = colores[2],
    by="runningScore",
    title = "Structural constituent of chromatin",
    pvalue_table = TRUE,
    base_size = 12,
    color.line = colores[1],
    color.vline = colores[5]
    ) +
    theme_minimal()
```

```{r,fig.width=5, fig.height=3}
# Buscar el ID real del término dentro del objeto GSEA
id <- enriquecimientos$Reactome@result |>
  filter(Description == "DNA methylation") |>
  pull(ID)

gseaplot(
    enriquecimientos$Reactome,
    geneSetID = id,
    color = colores[2],
    by="runningScore",
    title = "DNA methylation",
    pvalue_table = TRUE,
    base_size = 12,
    color.line = colores[1],
    color.vline = colores[5]
    ) +
    theme_minimal()
```

```{r,fig.width=5, fig.height=3}
# Buscar el ID real del término dentro del objeto GSEA
id <- enriquecimientos$Reactome@result |>
  filter(Description == "Interleukin-1 signaling") |>
  pull(ID)

gseaplot(
    enriquecimientos$Reactome,
    geneSetID = id,
    color = colores[2],
    by="runningScore",
    title = "Interleukin-1 signaling",
    pvalue_table = TRUE,
    base_size = 12,
    color.line = colores[1],
    color.vline = colores[5]
    ) +
    theme_minimal()
```

```{r,fig.width=5, fig.height=3}
# Buscar el ID real del término dentro del objeto GSEA
id <- enriquecimientos$Reactome@result |>
  filter(Description == "Interferon alpha/beta signaling") |>
  pull(ID)

gseaplot(
    enriquecimientos$Reactome,
    geneSetID = id,
    color = colores[2],
    by="runningScore",
    title = "Interferon alpha/beta signaling",
    pvalue_table = TRUE,
    base_size = 12,
    color.line = colores[1],
    color.vline = colores[5]
    ) +
    theme_minimal()
```

```{r,fig.width=5, fig.height=3}
# Buscar el ID real del término dentro del objeto GSEA
id <- enriquecimientos$GO@result |>
  filter(Description == "fibroblast proliferation") |>
  pull(ID)

gseaplot(
    enriquecimientos$GO,
    geneSetID = id,
    color = colores[2],
    by="runningScore",
    title = "fibroblast proliferation",
    pvalue_table = TRUE,
    base_size = 12,
    color.line = colores[1],
    color.vline = colores[5]
    ) +
    theme_minimal()
```


#### Tabla enriquecimientos genes

```{r}
# plot_df_gsea |> 
#     dplyr::select(Description, NES, pvalue, p.adjust, db) |>
#     write.csv("../../../../Artículos/HERVs/version0_no20E_RINadj/tablas/SupplementaryTable2.csv")
```


ddsSalmon_filtered_GC <- readRDS("../data/dds_DEG_filtered_GC.rds")
ddsSalmon_filtered_GC <- DESeq(ddsSalmon_filtered_GC)
res_DEG_uni_GC <- results(ddsSalmon_filtered_GC, contrast = c("condition", "COVID", "Control"))

## Annotation of results table
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

genes_DE_anotados |>
    filter(abs(log2FoldChange) >=1 & padj <= 0.05) |>
    dplyr::select(SYMBOL, log2FoldChange, padj) |> write.csv("../../../../Artículos/HERVs/tablas/SupplementaryTable1.csv")

## Housekeeping

hk_candidates <- intersect(
  genes_DE_anotados$SYMBOL,
  c(
      "ACTB", "APRT", "EEF1A1", "UBC", "HPRT1", "SDHA", "TUBB", "PPIA","TBP","PGK1","GUSB","PSMB4","VCP","POLR2A", "B2M", "RPLP1", "RPLP2", "TFRC", "ALAS1", "PUM1",
        "YWHAZ", "IPO8",
  "UBE2D2", "C1orf43", "EMC7"
  )
)

Genes_mcount_vst <- vst(assay(ddsSalmon_filtered_GC))
genes_entrezid_hk <- genes_DE_anotados |>
    filter(SYMBOL %in% hk_candidates) |> pull(ENTREZID)

mat <- Genes_mcount_vst[rownames(Genes_mcount_vst) %in% genes_entrezid_hk,]

df <- mat |>
    t() |> 
    as.data.frame() |> 
    rownames_to_column("sample") |> 
    mutate(across(where(is.factor), as.character)) |> 
    pivot_longer(cols = -sample) |>
    dplyr::left_join(patients_data_all_cov_RIN, by = c("sample" = "expediente")) |>
    dplyr::left_join(genes_DE_anotados, by = c("name" = "ENTREZID"))

## Functional enrichment

### gsea

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
    slice_max(
        order_by = abs(stat),
        n = 1,
        with_ties = FALSE
    ) |>
    ungroup()

geneList <- lista_GSEA.df$stat
names(geneList) <- lista_GSEA.df$ENTREZID

geneList <- sort(geneList, decreasing = TRUE)

set.seed(123)

enriquecimientos_GSEA <- list()


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

enriquecimientos_filtered_df <- data.frame()
enriquecimientos_filtrados <- enriquecimientos

for (i in names(enriquecimientos)) {
    
    df_filtered <- enriquecimientos[[i]]@result 
    
    enriquecimientos_filtered_df <- df_filtered |> 
        dplyr::select(Description, NES, p.adjust, pvalue) |>
        mutate(db = i) |> 
        rbind(enriquecimientos_filtered_df)
    
    
}

enriquecimientos_filtered_df <- enriquecimientos_filtered_df |> filter(p.adjust <= 0.05)

plot_df_gsea |> 
     dplyr::select(Description, NES, pvalue, p.adjust, db) |>
     write.csv("../../../../Artículos/HERVs/tablas/SupplementaryTable2.csv")

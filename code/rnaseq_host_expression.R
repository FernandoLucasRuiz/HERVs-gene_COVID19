patients_data_all <- readRDS("patients_data_all.rds")

temp <- read.delim("data/Salmon_results_gcBias/COV1C_quant_gc/quant.sf")
Tx2Gene <- AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg38.knownGene, keys = as.vector(temp[, 1]),
                                 keytype = "TXNAME", columns = c("GENEID", "TXNAME"))
Tx2Gene <- Tx2Gene[!is.na(Tx2Gene$GENEID), ]


salmonQ_GC <- dir("data/Salmon_results_gcBias/", 
                  recursive = T, 
                  pattern = "quant.sf", 
                  full.names = T)

names(salmonQ_GC) <- sub("_quant_gc/quant.sf","", sub("data/Salmon_results_gcBias//","", salmonQ_GC))

salmonCounts_GC <- tximport(salmonQ_GC, 
                            type = "salmon", 
                            tx2gene = Tx2Gene)

colnames(salmonCounts_GC$abundance) <- gsub("/", "", colnames(salmonCounts_GC$abundance))
colnames(salmonCounts_GC$counts) <- gsub("/", "", colnames(salmonCounts_GC$counts))
patients_data_Salmon_GC <- patients_data_all[colnames(salmonCounts_GC$abundance),]

ddsSalmon_import_GC  <- DESeqDataSetFromTximport(salmonCounts_GC, 
                                                 colData = patients_data_Salmon_GC, 
                                                 design = ~ condition)

keep_v02_S_GC <- rowSums(counts(ddsSalmon_import_GC)>10) >= 19
ddsSalmon_filtered_GC <- ddsSalmon_import_GC[keep_v02_S_GC,]

ddsSalmon_filtered_GC <- DESeq(ddsSalmon_filtered_GC)
res_DEG_uni_GC <- results(ddsSalmon_filtered_GC)

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


genes_DE_annotated <- res_DEG_uni_ann_GC

rm(res_DEG_uni_ann_GC)

## Housekeeping

hk_candidates <- intersect(
    genes_DE_annotated$SYMBOL,
    c(
        "ACTB", "APRT", "EEF1A1", "UBC", "HPRT1", "SDHA", "TUBB", "PPIA","TBP","PGK1","GUSB","PSMB4","VCP","POLR2A", "B2M", "RPLP1", "RPLP2", "TFRC", "ALAS1", "PUM1"
    )
    )

## heatmap
vsd <- vst(ddsSalmon_filtered_GC, blind = FALSE)
mat <- assay(vsd)

df <- as.data.frame(mat) |>
    rownames_to_column("genes") |>
    pivot_longer(-genes, names_to = "sample", values_to = "expression")

metadata <- as.data.frame(colData(ddsSalmon_filtered_GC)) |>
    rownames_to_column("sample")

df <- df |>
    left_join(metadata, by = "sample") 

tidyheatmap(
    df,
    rows = genes,
    columns = sample,
    values = expression,
    scale = "row", 
    color_legend_min = -1,
    color_legend_max = 1,
    colors = c("cadetblue3","#ffffff","#A52A2A"),
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

## umap
vsd <- vst(ddsSalmon_filtered_GC, blind = FALSE)
mat <- t(assay(vsd))

set.seed(1234)
um <- uwot::umap(mat)

## Functional enrichment

enrichments <- list()

lista_enriquecer.df <- genes_DE_annotated |>
    filter(abs(log2FoldChange) >= 1.5 & padj <= 0.05)
    
ego <- enrichGO(gene = unique(lista_enriquecer.df$ENTREZID),
                        keyType       = "ENTREZID",
                        OrgDb         = org.Hs.eg.db,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.1,
                        readable      = TRUE)
    
enrichments[["GO"]] <- ego
    
kk <- enrichKEGG(gene         = unique(lista_enriquecer.df$ENTREZID), 
                     organism     = 'hsa',
                     pvalueCutoff = 0.1)

enrichments[["KEGG"]] <- kk

wp <- enrichWP(gene = unique(lista_enriquecer.df$ENTREZID),
         organism = "Homo sapiens") 
    
enrichments[["WikiPathways"]] <- wp

x <- enrichPathway(gene=unique(lista_enriquecer.df$ENTREZID), 
                   pvalueCutoff = 0.1, 
                   readable=TRUE)

enrichments[["Reactome"]] <- x

x <- enrichDO(gene          = unique(lista_enriquecer.df$ENTREZID),
              ont           = "HDO",
              pvalueCutoff  = 0.1,
              pAdjustMethod = "BH",
              universe      = names(lista_enriquecer.df$ENTREZID),
              minGSSize     = 10,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)

enrichments[["DOSE_Disease"]] <- x

x <- enrichDO(gene          = unique(lista_enriquecer.df$ENTREZID),
              ont           = "HPO",
              pvalueCutoff  = 0.1,
              pAdjustMethod = "BH",
              universe      = names(lista_enriquecer.df$ENTREZID),
              minGSSize     = 10,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)

enrichments[["DOSE_Phenotype"]] <- x

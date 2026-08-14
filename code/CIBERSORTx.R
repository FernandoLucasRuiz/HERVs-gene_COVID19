# Preparing signature matrix
sce <- readH5AD(
    "../HLCA/dbb5ad81-1713-4aee-8257-396fbabe7c6e.h5ad",
    use_hdf5 = TRUE
)

metadata <- colData(sce) |>
    as.data.frame() 

metadata <- metadata |> 
    dplyr::filter(tissue %in% c("lung")) |>
    dplyr::filter(disease %in% c("normal", "COVID-19"))

sce <- sce[, colnames(sce) %in% rownames(metadata)]

obj <- as.Seurat(
    sce,
    counts = "soupX",
    data = NULL
)

metadata <- colData(sce) |>
    as.data.frame() 


metadata <- metadata |> 
    dplyr::filter(tissue %in% c("lung")) |>
    dplyr::filter(disease %in% c("normal", "COVID-19"))

sce <- sce[, colnames(sce) %in% rownames(metadata)]

set.seed(42)

cells_to_keep <- colData(sce) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_barcode") %>%
    filter(ann_level_2 != "Unknown") |>
    filter(ann_level_2 != "Hematopoietic stem cells") |>
    filter(ann_level_2 != "Submucosal Gland") |>
    filter(ann_level_2 != "Mesothelium") |>
  group_by(disease, ann_level_2) %>%
  slice_sample(prop = 1) %>%
  slice_head(n = 500) %>% 
  pull(cell_barcode)

sce_filtered <- sce[, cells_to_keep]

counts_matrix <- assay(sce_filtered, "X")

celltype_vector <- sce_filtered$ann_level_2 
names(celltype_vector) <- colnames(sce_filtered) 
disease <- sce_filtered$disease

label_df <- data.frame(
    Cell = names(celltype_vector),
    Label = factor(paste0(celltype_vector, "_", disease))
)

rownames(label_df) <- NULL

all(colnames(counts_matrix) %in% label_df$Cell)
all(label_df$Cell %in% colnames(counts_matrix))

scdata_sub <- as.matrix(counts_matrix)

idx <- match(colnames(scdata_sub), label_df$Cell)

scdata_renamed <- scdata_sub
colnames(scdata_renamed) <- label_df$Label[idx]

mat_out <- cbind(GeneSymbol = rownames(scdata_renamed), scdata_renamed)

ens_ids <- mat_out[, "GeneSymbol"]

ens_ids_clean <- sub("\\..*$", "", ens_ids)

symbols <- mapIds(
    org.Hs.eg.db,
    keys = ens_ids_clean,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
    )

mat_out[, "GeneSymbol"] <- symbols

mat_out_symbol <- mat_out[!is.na(mat_out[, "GeneSymbol"]), ]


write.table(mat_out_symbol, file = "scdata_sub_lung_HLCA_ann2_covid_normal.txt",
            sep = "\t", row.names = F, quote = FALSE, col.names = T)

# DEGs immune healthy vs COVID

pb_counts <- read.csv(
    "../HLCA/HLCA_immune_pseudobulk_counts_by_donor.csv.gz",
    row.names = 1,
    check.names = FALSE
    ) |>
    as.matrix()

pb_meta <- read.csv(
    "../HLCA/HLCA_immune_pseudobulk_metadata_by_donor.csv",
    row.names = 1,
    check.names = FALSE
    )

pb_meta$sample_id_pb <- rownames(pb_meta)

pb_meta <- pb_meta[
  match(colnames(pb_counts), pb_meta$sample_id_pb),
  ,
  drop = FALSE
]

storage.mode(pb_counts) <- "integer"

dds_pb <- DESeqDataSetFromMatrix(
    countData = pb_counts,
    colData   = pb_meta,
    design    = ~ condition
    )

keep <- rowSums(counts(dds_pb)) >= 10

dds_pb_filtered <- dds_pb[keep, ]

dds_pb_filtered <- DESeq(
    dds_pb_filtered,
    sfType = "poscounts"
    )

res_pb <- results(
    dds_pb_filtered,
    contrast = c("condition", "COVID", "Control"),
    alpha = 0.05
    )

res_pb_df <- as.data.frame(res_pb) |>
  rownames_to_column("gene_id") |>
  mutate(
    gene_id_clean = sub("\\..*$", "", gene_id),

    SYMBOL = mapIds(
      org.Hs.eg.db,
      keys = gene_id_clean,
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )[gene_id_clean],
    
    ENTREZID = mapIds(
      org.Hs.eg.db,
      keys = gene_id_clean,
      keytype = "ENSEMBL",
      column = "ENTREZID",
      multiVals = "first"
    )[gene_id_clean]
  ) |>
  relocate(SYMBOL, .after = gene_id)

res_bulk_df <- genes_DE_anotados |>
  dplyr::select(
    SYMBOL,
    bulk_log2FC = log2FoldChange,
    bulk_pvalue = pvalue,
    bulk_padj   = padj,
    bulk_stat   = stat
  )

res_pb_df <- res_pb_df |> drop_na() |>
    rename("logFC" ="log2FoldChange",
           "adj.P.Val"= "padj")

top_covid <- res_pb_df |>
  filter(adj.P.Val < 0.05, logFC > 1.5) |>
  arrange(adj.P.Val) |> head(100)

lista_enriquecer.df <- top_covid |>
  pull(ENTREZID) |>
  na.omit() |>
  unique()

enriquecimientos <- list()
    
ego <- enrichGO(gene = lista_enriquecer.df,
                        keyType       = "ENTREZID",
                        OrgDb         = org.Hs.eg.db,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.1,
                        readable      = TRUE)
    
enriquecimientos[["GO"]] <- ego
    
kk <- enrichKEGG(gene         = lista_enriquecer.df, 
                     organism     = 'hsa',
                     pvalueCutoff = 0.1)

enriquecimientos[["KEGG"]] <- kk

wp <- enrichWP(gene = lista_enriquecer.df,
         organism = "Homo sapiens") 
    
enriquecimientos[["WikiPathways"]] <- wp

x <- enrichPathway(gene=lista_enriquecer.df,
                   pvalueCutoff = 0.1, 
                   readable=TRUE)

enriquecimientos[["Reactome"]] <- x

enriquecimientos_filtered_df <- data.frame()
enriquecimientos_filtrados <- enriquecimientos_pseudobulk

for (i in names(enriquecimientos)) {
    
    df_filtered <- enriquecimientos[[i]]@result 
    
    enriquecimientos_filtered_df <- df_filtered |> 
        dplyr::select(Description, Count, GeneRatio, p.adjust, pvalue, geneID) |>
        mutate(db = i) |> 
        rbind(enriquecimientos_filtered_df)
    
    
}

enriquecimientos_filtered_df_pseudobulk <- enriquecimientos_filtered_df |> filter(p.adjust <= 0.05)

df_enriquecido |>
    dplyr::select(Description, GeneRatio, pvalue, p.adjust) |>
    write.csv("../../../../Artículos/HERVs/tablas/SupplementaryTable3.csv")

genes_pseudobulk <- res_pb_df |>
    filter(adj.P.Val < 0.05) |>
    filter(logFC > -10) |>
    filter(SYMBOL != "SMG1P7") |>
    arrange(adj.P.Val) |>
    drop_na() |>
    dplyr::select(SYMBOL, logFC, adj.P.Val) |>
    rename("logFC_pseudobulk" = "logFC" ,
           "padj_pseudobulk" = "adj.P.Val" )

comparison_pb_bulk <- genes_DE_anotados |>
    filter(padj < 0.05) |>
    dplyr::select(SYMBOL, log2FoldChange, padj) |>
    rename("logFC_bulk" = "log2FoldChange" ,
           "padj_bulk" = "padj" ) |>
    inner_join(genes_pseudobulk, by= "SYMBOL")  |> 
    filter(
        is.finite(logFC_bulk),
        is.finite(logFC_pseudobulk)
    )


spearman_test <- cor.test(
    comparison_pb_bulk$logFC_bulk,
    comparison_pb_bulk$logFC_pseudobulk,
    method = "spearman",
    exact = FALSE
    )

# DEGs Fibroblast healthy vs COVID

pb_counts <- read.csv(
  "../HLCA/HLCA_Fibroblast_pseudobulk_counts_by_donor.csv.gz",
  row.names = 1,
  check.names = FALSE
) |>
  as.matrix()

pb_meta <- read.csv(
  "../HLCA/HLCA_Fibroblast_pseudobulk_metadata_by_donor.csv",
  row.names = 1,
  check.names = FALSE
)

pb_meta$sample_id_pb <- rownames(pb_meta)

# Ordenar los metadatos igual que las columnas de la matriz
pb_meta <- pb_meta[
  match(colnames(pb_counts), pb_meta$sample_id_pb),
  ,
  drop = FALSE
]

storage.mode(pb_counts) <- "integer"

dds_pb <- DESeqDataSetFromMatrix(
  countData = pb_counts,
  colData   = pb_meta,
  design    = ~ condition
)

keep <- rowSums(counts(dds_pb)) >= 10

dds_pb_filtered <- dds_pb[keep, ]


dds_pb_filtered <- DESeq(
    dds_pb_filtered,
    sfType = "poscounts"
    )

res_pb <- results(
    dds_pb_filtered,
    contrast = c("condition", "COVID", "Control"),
    alpha = 0.05
    )

res_pb_df <- as.data.frame(res_pb) |>
  rownames_to_column("gene_id") |>
  mutate(
    gene_id_clean = sub("\\..*$", "", gene_id),

    SYMBOL = mapIds(
      org.Hs.eg.db,
      keys = gene_id_clean,
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )[gene_id_clean],
    
    ENTREZID = mapIds(
      org.Hs.eg.db,
      keys = gene_id_clean,
      keytype = "ENSEMBL",
      column = "ENTREZID",
      multiVals = "first"
    )[gene_id_clean]
  ) |>
  relocate(SYMBOL, .after = gene_id)

res_bulk_df <- genes_DE_anotados |>
  dplyr::select(
    SYMBOL,
    bulk_log2FC = log2FoldChange,
    bulk_pvalue = pvalue,
    bulk_padj   = padj,
    bulk_stat   = stat
  )

res_pb_df <- res_pb_df |> drop_na() |>
    rename("logFC" ="log2FoldChange",
           "adj.P.Val"= "padj")

top_covid <- res_pb_df |>
  filter(adj.P.Val < 0.05, logFC > 1.5) |>
  arrange(adj.P.Val) |> head(100)

lista_enriquecer.df <- top_covid |>
  pull(ENTREZID) |>
  na.omit() |>
  unique()

enriquecimientos <- list()
    
ego <- enrichGO(gene = lista_enriquecer.df,
                        keyType       = "ENTREZID",
                        OrgDb         = org.Hs.eg.db,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.1,
                        readable      = TRUE)
    
enriquecimientos[["GO"]] <- ego
    
kk <- enrichKEGG(gene         = lista_enriquecer.df, 
                     organism     = 'hsa',
                     pvalueCutoff = 0.1)

enriquecimientos[["KEGG"]] <- kk

wp <- enrichWP(gene = lista_enriquecer.df,
         organism = "Homo sapiens") 
    
enriquecimientos[["WikiPathways"]] <- wp

x <- enrichPathway(gene=lista_enriquecer.df,
                   pvalueCutoff = 0.1, 
                   readable=TRUE)

enriquecimientos[["Reactome"]] <- x

enriquecimientos_filtered_df <- data.frame()
enriquecimientos_filtrados <- enriquecimientos

for (i in names(enriquecimientos)) {
    
    df_filtered <- enriquecimientos[[i]]@result 
    
    enriquecimientos_filtered_df <- df_filtered |> 
        dplyr::select(Description, Count, GeneRatio, p.adjust, pvalue, geneID) |>
        mutate(db = i) |> 
        rbind(enriquecimientos_filtered_df)
    
    
}

enriquecimientos_filtered_df_pseudobulk <- enriquecimientos_filtered_df |> filter(p.adjust <= 0.05)

genes_pseudobulk <- res_pb_df |>
    filter(adj.P.Val < 0.05) |>
    filter(logFC > -10) |>
    arrange(adj.P.Val) |>
    drop_na() |>
    dplyr::select(SYMBOL, logFC, adj.P.Val) |>
    rename("logFC_pseudobulk" = "logFC" ,
           "padj_pseudobulk" = "adj.P.Val" )

comparison_pb_bulk <- genes_DE_anotados |>
    filter(padj < 0.05) |>
    dplyr::select(SYMBOL, log2FoldChange, padj) |>
    rename("logFC_bulk" = "log2FoldChange" ,
           "padj_bulk" = "padj" ) |>
    inner_join(genes_pseudobulk, by= "SYMBOL")  |> 
    filter(
        is.finite(logFC_bulk),
        is.finite(logFC_pseudobulk)
    )

spearman_test <- cor.test(
    comparison_pb_bulk$logFC_bulk,
    comparison_pb_bulk$logFC_pseudobulk,
    method = "spearman",
    exact = FALSE
    )

# Cibersorting

mat_genes_todos <- counts(ddsSalmon_filtered_GC, normalized = TRUE)


mat_genes_todos |> as.data.frame() |> rownames_to_column("ENTREZID") |>
    left_join(genes_DE_anotados |> dplyr::select(ENTREZID, SYMBOL), by="ENTREZID") |>
    filter(!is.na(SYMBOL)) |>
    column_to_rownames("SYMBOL") |> dplyr::select(c(-ENTREZID)) |> rownames_to_column("gene_name") |>
    write.table(file = "todos_genes.txt",
            sep = "\t", row.names = F, quote = FALSE, col.names = T)


matriz_fractions <- read.table("../HLCA/docker_results/todos_ann2_covid-normal_results//CIBERSORTx_Adjusted.txt", sep = "\t", header = T, stringsAsFactors = FALSE)

cibersort_fraction_2 <- matriz_fractions |> 
    dplyr::select(-c("P.value", "RMSE", Correlation))|>
    as.data.frame() |> 
    rename(expediente = "Mixture") 

df <- inner_join(patients_data_all_cov_RIN, cibersort_fraction_2, by = "expediente")
df_plot <- df |>
    dplyr::select(condition, colnames(cibersort_fraction_2)[-1]) 


patients_data <- readRDS("patients_data.rds")

files <- list.files(path = "data/HERVs/Telescope_may2024/", 
             pattern = "-TE_counts.tsv", 
             all.files = FALSE,full.names = FALSE)

file_names <- sub("-TE_counts.tsv","", files)
matrixcount_all <- data.frame(transcript=character())


for (i in seq(length(file_names))) {
  file_bam <- read.delim(paste("data/HERVs/Telescope_may2024/",files[i], sep ="")) 
  names(file_bam)[2] <-  file_names[i] 
  matrixcount_all <- merge(matrixcount_all,file_bam,by="transcript",all=TRUE) 
}


rownames(matrixcount_all) <- matrixcount_all$transcript

colnames(matrixcount_all) <- sub("_retro","",colnames(matrixcount_all))

matrixcount_all <- matrixcount_all[rownames(matrixcount_all)!="__no_feature",]

patients_data_all <- read.table("data/HERVs/groups_patients_all.txt")
patients_data_all$condition <- factor(patients_data_all$condition, levels=c("Control", "COVID"))

matrixcount_complete <- matrixcount_all[,rownames(patients_data_all)]

matrixcount_complete[is.na(matrixcount_complete)] <- 0

dds_all <- DESeqDataSetFromMatrix(countData = matrixcount_complete,
                              colData = patients_data,
                              design = ~ condition)


min_counts <- 10
n_samples <- min(table(patients_data$condition))
keep_v02_A <- rowSums(counts(dds_all)>min_counts) >= n_samples
dds_prefil_all <- dds_all[keep_v02_A,]

dds_DE_all <- DESeq(dds_prefil_all) #DE analysis
res_HERVs_all <- results(dds_DE_all) #results

HERVs_DE_all <- res_HERVs_all |>
    as.data.frame() |>
    rownames_to_column("HERVs") |>
    mutate(
        expression= case_when(
            log2FoldChange < 0 & padj <= 0.05 ~ "Significative downregulated",
            log2FoldChange > 0 & padj < 0.05 ~ "Significative upregulated",
            T ~ "No significative"
        )
    )

rowranges <- read.delim("data/TE_annotation.v2.0.tsv", 
                        sep = "\t", 
                        header = TRUE, 
                        quote = "", 
                        row.names = 1)

rowranges$locus <- rownames(rowranges)

HERVs_family_counts <- read.delim("data/HERVs_family_counts.txt", header=FALSE)

names(HERVs_family_counts) <- c("Family", "counts", "Family_group")

HERVs_annotation <- merge(rowranges, HERVs_family_counts[,c("Family", "Family_group")], all=TRUE)

rownames(HERVs_annotation) <- HERVs_annotation$locus

HERVs_annotation <- HERVs_annotation |> 
    rownames_to_column("HERVs")

HERVs_DE_ann_all <- merge(as.data.frame(HERVs_DE_all), HERVs_annotation, by='HERVs')
rownames(HERVs_DE_ann_all) <- HERVs_DE_ann_all$Row.names

### Hetamap
vsd <- vst(dds_DE_all, blind = FALSE)
mat <- assay(vsd)

df <- as.data.frame(mat) |>
    rownames_to_column("HERVs") |>
    pivot_longer(-HERVs, names_to = "sample", values_to = "expression")

metadata <- as.data.frame(colData(dds_DE_all)) |>
    rownames_to_column("sample")

df <- df |>
    left_join(metadata, by = "sample") |>
    left_join(HERVs_DE_annotated, by = "HERVs") |>
    filter(Class == "HERV")

tidyheatmap(
    df,
    rows = HERVs,
    columns = sample,
    values = expression.x,
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
        condition = c(Control= "darkolivegreen4", COVID ="#FFA500")
        )
    ) 

### UMAP
vsd <- vst(dds_DE_all, blind = FALSE)
mat <- t(assay(vsd))

mat <- mat[,colnames(mat) %in% HERVs_DE_annotated$HERVs]

set.seed(1234)
um <- uwot::umap(mat)

## functional enrichment gene associated with Hervs

df <- HERVs_DE_annotated |>
    filter(padj <=0.05 & log2FoldChange > 1.5) |>
    dplyr::select(HERVs, IntersectedGeneID, ClosestDownstream_id, ClosestDownstream_id) |>
    pivot_longer(-HERVs) |>
    drop_na() |> filter(value != "None") |>
    separate(value, into = "value", sep = ",")

genes_asociados_HERVs_upregulados <- df |>
    pull(value)

Gene_ann_GC <- AnnotationDbi::select(org.Hs.eg.db,
                 keys = genes_asociados_HERVs_upregulados,
                 keytype = "ENSEMBL",
                 columns = c("SYMBOL", "ENTREZID"))

res_DEG_uni_ann_GC <- inner_join(df, Gene_ann_GC, by =c("value" = "ENSEMBL"))

enriquecimientos <- list()
    
ego <- enrichGO(gene = unique(res_DEG_uni_ann_GC$ENTREZID),
                        keyType       = "ENTREZID",
                        OrgDb         = org.Hs.eg.db,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.1,
                        readable      = TRUE)
    
enriquecimientos[["GO"]] <- ego
    
kk <- enrichKEGG(gene       = res_DEG_uni_ann_GC$ENTREZID, 
                     organism     = 'hsa',
                     pvalueCutoff = 0.1)

enriquecimientos[["KEGG"]] <- kk

wp <- enrichWP(gene = unique(res_DEG_uni_ann_GC$ENTREZID),
         organism = "Homo sapiens") 
    
enriquecimientos[["WikiPathways"]] <- wp

x <- enrichPathway(gene=unique(res_DEG_uni_ann_GC$ENTREZID), 
                   pvalueCutoff = 0.1, 
                   readable=TRUE)

enriquecimientos[["Reactome"]] <- x

x <- enrichDO(gene          = unique(res_DEG_uni_ann_GC$ENTREZID),
              ont           = "HDO",
              pvalueCutoff  = 0.1,
              pAdjustMethod = "BH",
              universe      = names(res_DEG_uni_ann_GC$ENTREZID),
              minGSSize     = 10,
              maxGSSize     = 500,
              qvalueCutoff  = 0.1,
              readable      = FALSE)

enriquecimientos[["DOSE_Disease"]] <- x

x <- enrichDO(gene          = unique(res_DEG_uni_ann_GC$ENTREZID),
              ont           = "HPO",
              pvalueCutoff  = 0.1,
              pAdjustMethod = "BH",
              universe      = names(res_DEG_uni_ann_GC$ENTREZID),
              minGSSize     = 10,
              maxGSSize     = 500,
              qvalueCutoff  = 0.1,
              readable      = FALSE)

enriquecimientos[["DOSE_Phenotype"]] <- x

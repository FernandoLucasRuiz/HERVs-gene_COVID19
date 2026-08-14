# HERVs differential expression analysis of COVID lung samples vs. health lung samples 

## Pre-filtering

dds_all <- DESeqDataSetFromMatrix(countData = matrixcount_complete,
                              colData = patients_data_all_cov_RIN,
                              design = ~ RIN_centered + condition)


min_counts <- 10
n_samples <- min(table(patients_data_all_cov_RIN$condition))
keep_v02_A <- rowSums(counts(dds_all)>min_counts) >= n_samples
dds_prefil_all <- dds_all[keep_v02_A,]

dds_DE_all <- DESeq(dds_prefil_all) #DE analysis
res_HERVs_all <- results(dds_DE_all, contrast = c("condition", "COVID", "Control")) #results

HERVs_DE_all <- res_HERVs_all |>
    as.data.frame() |>
    rownames_to_column("HERVs") 

rowranges <- read.delim("../data/TE_annotation.v2.0.tsv", 
                        sep = "\t", 
                        header = TRUE, 
                        quote = "", 
                        row.names = 1)

rowranges$locus <- rownames(rowranges)

HERVs_family_counts <- read.delim("../data/HERVs_family_counts.txt", header=FALSE)

names(HERVs_family_counts) <- c("Family", "counts", "Family_group")


HERVs_annotation <- merge(rowranges, HERVs_family_counts[,c("Family", "Family_group")], all=TRUE)

rownames(HERVs_annotation) <- HERVs_annotation$locus

HERVs_annotation <- HERVs_annotation |> 
    rownames_to_column("HERVs")

HERVs_DE_ann_all <- merge(as.data.frame(HERVs_DE_all), HERVs_annotation, by='HERVs')
rownames(HERVs_DE_ann_all) <- HERVs_DE_ann_all$Row.names

HERVs_DE_anotados <- HERVs_DE_ann_all |>
    filter(Class == "HERV")

rm(HERVs_DE_ann_all)

df <- HERVs_DE_anotados |>
    filter(padj <=0.05 & log2FoldChange > 1.5) |>
    dplyr::select(HERVs, IntersectedGeneID, ClosestDownstream_id, ClosestUpstream_id) |>
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


enriquecimientos_filtered_df <- data.frame()

for (i in names(enriquecimientos)) {
    
    df_filtered <- enriquecimientos[[i]]@result 
    
    enriquecimientos_filtered_df <- df_filtered |> 
        dplyr::select(Description, GeneRatio, Count, pvalue, p.adjust, geneID) |>
        mutate(db = i) |> 
        rbind(enriquecimientos_filtered_df)

}

enriquecimientos_filtered_df <- enriquecimientos_filtered_df |>
    filter(pvalue <= 0.05)

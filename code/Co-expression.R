## Normalizing matrix counts con DEG and DHERVs

dds_prefil_all <-readRDS("../data/no_20E_RINadj/dds_HERVs_filtered.rds")
dds_DE_all <- DESeq(dds_prefil_all) #DE analysis
res_all <- results(dds_DE_all, contrast = c("condition", "COVID", "Control")) #results

ddsSalmon_filtered_GC <- readRDS("../data/no_20E_RINadj/dds_DEG_filtered_GC.rds")
ddsSalmon_filtered_GC <- DESeq(ddsSalmon_filtered_GC)
res_DEG_uni_GC <- results(ddsSalmon_filtered_GC, contrast = c("condition", "COVID", "Control"))

HERVs_DE_ann_all <- readRDS("../data/no_20E_RINadj/HERVs_DE_ann_all.rds")
HERVs_DE_anotados <- HERVs_DE_ann_all |>
    filter(Class == "HERV")

rm(HERVs_DE_ann_all)


dds_herv <- dds_DE_all

herv_counts <- DESeq2::counts(dds_herv)

gene_counts <- DESeq2::counts(ddsSalmon_filtered_GC)

gene_counts <- gene_counts[
  ,
  colnames(herv_counts),
  drop = FALSE
]

min_samples <- min(table(colData(dds_herv)$condition))

# Filtrado por expresión
keep_herv <- rowSums(herv_counts >= 5) >= min_samples
keep_gene <- rowSums(gene_counts >= 5) >= min_samples

# VST
herv_vsd <- DESeq2::varianceStabilizingTransformation(
    dds_herv[keep_herv, ],
    blind = FALSE
    )

gene_vsd <- DESeq2::varianceStabilizingTransformation(
    ddsSalmon_filtered_GC[keep_gene, ],
    blind = FALSE
    )

HERVs_mcount_vst <- SummarizedExperiment::assay(herv_vsd)
Genes_mcount_vst <- SummarizedExperiment::assay(gene_vsd)

Genes_mcount_vst <- Genes_mcount_vst[
    ,
    colnames(HERVs_mcount_vst),
    drop = FALSE
    ]

Genes_HERVs_mcounts_vst <- rbind(
    HERVs_mcount_vst,
    Genes_mcount_vst
    )

Genes_HERVs_mcounts_vst_t <- as.data.frame(
    t(Genes_HERVs_mcounts_vst)
    )

Genes_HERVs_filtered <- GWENA::filter_low_var(
    Genes_HERVs_mcounts_vst_t,
    pct = 0.7,
    type = "median"
    )

## Network building

threads_to_use <- 4
net_DE <- build_net(
    Genes_HERVs_filtered,
                    n_threads = threads_to_use, 
                    fit_cut_off=0.8, 
                    cor_func = "spearman",
                    keep_matrices = "cor"
                    )

## Modules detection

modules <- detect_modules(
    Genes_HERVs_filtered,
    net_DE$network, 
    detailled_result = TRUE,
    merge_threshold = 0.25)

## Phenotypic association

patients_data_all_GWENA <- patients_data_all_cov_RIN
patients_data_all_GWENA$condition <- as.character(patients_data_all_GWENA$condition)

patients_data_all_GWENA_sex <- patients_data_all_cov_RIN[rownames(Genes_HERVs_filtered),c("expediente", "condition", "sexo")]

patients_data_all_GWENA_sex$condition <- as.character(patients_data_all_GWENA_sex$condition)
patients_data_all_GWENA_sex$sexo <- as.character(patients_data_all_GWENA_sex$sexo)

phenotype_association <- associate_phenotype(
    modules$modules_eigengenes, 
    (patients_data_all_GWENA_sex %>% dplyr::select(c(condition, sexo))))



df <- phenotype_association$association |> dplyr::select(COVID)
M <- as.matrix(df)

df2 <- phenotype_association$pval |> dplyr::select(COVID)
testRes <- as.matrix(df2)


get_kME_df <- function(module_genes, ME_col, module_name) {
  
  expr_module <- as.matrix(
    Genes_HERVs_filtered[, module_genes, drop = FALSE]
    
  )
  
  common_samples <- intersect(rownames(expr_module), rownames(ME))
  
  expr_module <- expr_module[common_samples, , drop = FALSE]
  ME_vector <- ME[common_samples, ME_col]
  
  kME <- cor(
    expr_module,
    ME_vector,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  
  kME_pvalue <- apply(
    expr_module,
    2,
    function(x) {
      cor.test(
        x,
        ME_vector,
        method = "spearman",
        exact = FALSE
      )$p.value
    }
  )
  
  data.frame(
    node = colnames(expr_module),
    kME = as.numeric(kME[, 1]),
    pvalue = as.numeric(kME_pvalue),
    padj = p.adjust(kME_pvalue, method = "BH"),
    modulo = module_name,
    stringsAsFactors = FALSE
  ) |>
    filter(
      node %in% module_genes,
      node %in% HERVs_DE_anotados$HERVs
    )
}

kME_df_M1 <- get_kME_df(
  module_genes = modules$modules$`1`,
  ME_col = "ME1",
  module_name = "ME1"
)

kME_df_M2 <- get_kME_df(
  module_genes = modules$modules$`2`,
  ME_col = "ME2",
  module_name = "ME2"
)

## M1

### Functional enrichment 

enriquecimientos <- list()

lista_enriquecer <- modules$modules$`1`[modules$modules$`1` %in% genes_DE_anotados$ENTREZID]
    
ego <- enrichGO(gene = unique(lista_enriquecer),
                        keyType       = "ENTREZID",
                        OrgDb         = org.Hs.eg.db,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.1,
                        readable      = TRUE)
    
enriquecimientos[["GO"]] <- ego
    
kk <- enrichKEGG(gene         = unique(lista_enriquecer), 
                     organism     = 'hsa',
                     pvalueCutoff = 0.1)

enriquecimientos[["KEGG"]] <- kk

wp <- enrichWP(gene = unique(lista_enriquecer),
         organism = "Homo sapiens") 
    
enriquecimientos[["WikiPathways"]] <- wp

x <- enrichPathway(gene=unique(lista_enriquecer), 
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

enriquecimientos_filtered_df <- enriquecimientos_filtered_df |> filter(p.adjust <= 0.05)

## M2

### Functional enrichment

enriquecimientos <- list()

lista_enriquecer <- modules$modules$`2`[modules$modules$`2` %in% genes_DE_anotados$ENTREZID]
    
ego <- enrichGO(gene = unique(lista_enriquecer),
                        keyType       = "ENTREZID",
                        OrgDb         = org.Hs.eg.db,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.1,
                        readable      = TRUE)
    
enriquecimientos[["GO"]] <- ego
    
kk <- enrichKEGG(gene         = unique(lista_enriquecer), 
                     organism     = 'hsa',
                     pvalueCutoff = 0.1)

enriquecimientos[["KEGG"]] <- kk

wp <- enrichWP(gene = unique(lista_enriquecer),
         organism = "Homo sapiens") 
    
enriquecimientos[["WikiPathways"]] <- wp

x <- enrichPathway(gene=unique(lista_enriquecer), 
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

enriquecimientos_filtered_df <- enriquecimientos_filtered_df |> filter(p.adjust <= 0.05)

### graph

module_example <- modules$modules$`2`
g <- build_graph_from_sq_mat(net_DE$network[module_example, module_example])

graph_nodes <- as.character(V(g)$name)

herv_pattern <- HERVs_DE_anotados$HERVs %>%
    as.character() %>%
    na.omit() %>%
    unique()

go_terms <- list(
    chromatin = c(
        "GO:0016568",  # Chromatin modification
        "GO:0006338",  # Chromatin remodeling
        "GO:0006306"   # DNA methylation
        )
    )

get_go_entrez <- function(go_ids) {

    AnnotationDbi::select(
        org.Hs.eg.db,
        keys = go_ids,
        keytype = "GOALL",
        columns = "ENTREZID"
        ) %>%
    pull(ENTREZID) %>%
    as.character() %>%
    na.omit() %>%
    unique()
    }


go_genes_entrez <- lapply(
    go_terms,
    get_go_entrez
    )


all_go_entrez <- go_genes_entrez %>%
    unlist(use.names = FALSE) %>%
    as.character() %>%
    unique()

go_id_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = all_go_entrez,
    keytype = "ENTREZID",
    columns = c(
        "SYMBOL",
        "ENSEMBL"
    )
    ) %>%
    mutate(
        ENTREZID = as.character(ENTREZID),
        SYMBOL = as.character(SYMBOL),
        ENSEMBL = as.character(ENSEMBL)
        ) %>%
    distinct(
        ENTREZID,
        SYMBOL,
        ENSEMBL
        )


go_genes <- lapply(
    go_genes_entrez,
    function(entrez_ids) {

        ids <- go_id_map %>%
            filter(
                ENTREZID %in% as.character(entrez_ids)
            ) %>%
        dplyr::select(
            ENTREZID,
            SYMBOL,
            ENSEMBL
            ) %>%
        unlist(use.names = FALSE) %>%
        as.character()

        ids <- ids[!is.na(ids) &ids != ""]

        unique(ids)
        }
    )



gene_to_symbol <- bind_rows(

    go_id_map %>%
        transmute(
            id = ENTREZID,
            label = SYMBOL
            ),

    go_id_map %>%
        transmute(
            id = SYMBOL,
            label = SYMBOL
            ),

    go_id_map %>%
        transmute(
            id = ENSEMBL,
            label = SYMBOL
            )
    ) %>%
    filter(
        !is.na(id),
        id != "",
        !is.na(label),
        label != ""
    ) %>%
    distinct(
        id,
        .keep_all = TRUE
    ) %>%
    tibble::deframe()



gene_id_to_entrez <- bind_rows(

    go_id_map %>%
        transmute(
            id = ENTREZID,
            ENTREZID = ENTREZID
            ),

    go_id_map %>%
        transmute(
            id = SYMBOL,
            ENTREZID = ENTREZID
            ),

    go_id_map %>%
        transmute(
            id = ENSEMBL,
            ENTREZID = ENTREZID
            )

    ) %>%
    filter(
        !is.na(id),
        id != "",
        !is.na(ENTREZID),
        ENTREZID != ""
        ) %>%
    distinct(
        id,
        .keep_all = TRUE
        ) %>%
    tibble::deframe()


gene_de <- genes_DE_anotados %>%
    transmute(
        ENTREZID = as.character(ENTREZID),
        GENE_log2FC = suppressWarnings(
        as.numeric(as.character(log2FoldChange))
        )) %>%
    filter(
        !is.na(ENTREZID),
        ENTREZID != "",
        !is.na(GENE_log2FC),
        GENE_log2FC > 0
        ) %>%
    arrange(
        desc(GENE_log2FC)
        ) %>%
    distinct(
        ENTREZID,
        .keep_all = TRUE
    )


edge_df <- igraph::as_data_frame(
    g,
    what = "edges"
    ) %>%
    mutate(
        from = as.character(from),
        to = as.character(to),
        weight = as.numeric(weight),
        abs_weight = abs(weight)
        )



hervs_in_graph <- intersect(graph_nodes,herv_pattern)

gene_nodes <- setdiff(graph_nodes,herv_pattern)


gene_node_keys <- sub("\\.[0-9]+$","",gene_nodes)


chromatin_genes_in_graph <- gene_nodes[gene_node_keys %in% go_genes$chromatin]


build_term_graph <- function(
        term_name,
        edge_df,
        go_genes,
        herv_pattern,
        gene_to_symbol,
        kme_df,
    gene_de,
    gene_id_to_entrez,
    top_hervs = 20,
    min_kme = 0.1,
    top_genes_per_herv = 5
) {

  term_gene_ids <- as.character(
    go_genes[[term_name]]
  )

  herv_gene_edges <- edge_df %>%
    mutate(
      from_is_herv = from %in% herv_pattern,
      to_is_herv = to %in% herv_pattern
    ) %>%
    filter(
      xor(
        from_is_herv,
        to_is_herv
      )
    ) %>%
    transmute(
      HERV = if_else(
        from_is_herv,
        from,
        to
      ),
      GENE = if_else(
        from_is_herv,
        to,
        from
      ),
      weight = weight,
      abs_weight = abs(weight)
    ) %>%
    distinct(
      HERV,
      GENE,
      .keep_all = TRUE
    )


  herv_gene_edges <- herv_gene_edges %>%
    mutate(
      GENE_key = sub(
        "\\.[0-9]+$",
        "",
        GENE
      )
    )


  term_edges <- herv_gene_edges %>%
    filter(
      GENE_key %in% term_gene_ids
    )


  term_edges <- term_edges %>%
    mutate(
      ENTREZID = unname(
        gene_id_to_entrez[GENE_key]
      )
    ) %>%
    left_join(
      gene_de,
      by = "ENTREZID"
    ) %>%
    filter(
      !is.na(GENE_log2FC),
      GENE_log2FC > 0
    ) %>%
    mutate(
      GENE_DE_direction = "Upregulated in COVID"
    )

  kme_hervs <- kme_df %>%
    transmute(
      HERV = as.character(node),
      kME = suppressWarnings(
        as.numeric(as.character(kME))
      )
    ) %>%
    filter(
      !is.na(HERV),
      !is.na(kME)
    ) %>%
    distinct(
      HERV,
      .keep_all = TRUE
    )


  connected_hervs <- unique(
    term_edges$HERV
  )


  n_overlap <- sum(
    connected_hervs %in% kme_hervs$HERV
  )


  kme_hervs <- kme_hervs %>%
    filter(
      HERV %in% connected_hervs,
      kME > 0
    )

  selected_hervs <- kme_hervs %>%
    filter(
      kME >= min_kme
    ) %>%
    arrange(
      desc(kME)
    ) %>%
    slice_head(
      n = top_hervs
    )
  

  term_edges <- term_edges %>%
    filter(
      HERV %in% selected_hervs$HERV
    ) %>%
    left_join(
      selected_hervs,
      by = "HERV"
    )


  if (!is.null(top_genes_per_herv)) {

    term_edges_small <- term_edges %>%
      group_by(HERV) %>%
      slice_max(
        order_by = abs_weight,
        n = top_genes_per_herv,
        with_ties = FALSE
      ) %>%
      ungroup()

  } else {

    term_edges_small <- term_edges
  }


  node_names <- unique(
    c(
      term_edges_small$HERV,
      term_edges_small$GENE
    )
  )


  nodes_df <- tibble(
    name = node_names
  ) %>%
    mutate(
      type = if_else(
        name %in% term_edges_small$HERV,
        "HERV",
        "GENE"
      ),

      lookup_id = if_else(
        type == "GENE",
        sub(
          "\\.[0-9]+$",
          "",
          name
        ),
        name
      ),

      mapped_symbol = unname(
        gene_to_symbol[lookup_id]
      ),

      label = case_when(
        type == "HERV" ~ name,

        !is.na(mapped_symbol) &
          mapped_symbol != "" ~ mapped_symbol,

        TRUE ~ name
      )
    ) %>%
    left_join(
      selected_hervs %>%
        rename(
          name = HERV
        ),
      by = "name"
    ) %>%
    dplyr::select(
      name,
      type,
      label,
      kME
    )


  g_term <- graph_from_data_frame(
    term_edges_small %>%
      dplyr::select(
        from = HERV,
        to = GENE,
        weight,
        abs_weight,
        GENE_log2FC,
        GENE_DE_direction
      ),
    directed = FALSE,
    vertices = nodes_df
  )

  node_match <- match(
    V(g_term)$name,
    nodes_df$name
  )


  V(g_term)$type_custom <- nodes_df$type[
    node_match
  ]

  V(g_term)$label_custom <- nodes_df$label[
    node_match
  ]

  V(g_term)$kME <- nodes_df$kME[
    node_match
  ]

  V(g_term)$degree_custom <- degree(
    g_term
  )


  E(g_term)$abs_weight <- abs(
    E(g_term)$weight
  )


  V(g_term)$strength_custom <- strength(
    g_term,
    weights = E(g_term)$abs_weight
  )


  return(g_term)
}


g_chromatin <- build_term_graph(
  term_name = "chromatin",
  edge_df = edge_df,
  go_genes = go_genes,
  herv_pattern = herv_pattern,
  gene_to_symbol = gene_to_symbol,
  kme_df = kME_df_M2,
  gene_de = gene_de,
  gene_id_to_entrez = gene_id_to_entrez,

  top_hervs = 3,
  min_kme = 0.9,

  top_genes_per_herv = 20
)

## Position in genome

mart_object_mygenes <- useEnsembl(biomart = "genes",dataset = "hsapiens_gene_ensembl")

my_genes = genes_DE_anotados$SYMBOL

mybdd <- getBM(mart=mart_object_mygenes, attributes=c('hgnc_symbol', 'description', 'chromosome_name','start_position', 'end_position', 'strand','ensembl_gene_id'),
            filters='hgnc_symbol', values=my_genes) # where df is a data.frame with all your requested info

### M1

DEG_chrom_FC <- merge(mybdd[,c("chromosome_name","start_position","end_position", "hgnc_symbol")], genes_DE_anotados, by.x="hgnc_symbol", by.y="SYMBOL")


DEG_chrom <- DEG_chrom_FC %>%
  filter(!grepl("^H", chromosome_name)) %>%
  mutate(chr=paste("chr",chromosome_name, sep="")) %>%
  dplyr::select("chr","start_position","end_position", "log2FoldChange","hgnc_symbol", "ENTREZID") %>%
  arrange(chr)


DEG_chrom_sort <- 
  data.frame(
  chr = DEG_chrom$chr,
  start = DEG_chrom$start_position,
  end = DEG_chrom$end_position,
  value = DEG_chrom$log2FoldChange,
  gene = DEG_chrom$hgnc_symbol,
  gene_entrez = DEG_chrom$ENTREZID
)


DEG_chrom_sort_up <-
  DEG_chrom_sort %>%
    filter(gene_entrez %in% c(modules$modules$`1`)) |>
    dplyr::select(-gene_entrez) |>
  filter(value>0)

DEG_chrom_sort_down <-
  DEG_chrom_sort %>%
    filter(gene_entrez %in% c(modules$modules$`1`)) |>
    dplyr::select(-gene_entrez) |>
  filter(value<0)

DEGs_chrom_list = list(DEG_chrom_sort_up, DEG_chrom_sort_down)

HERVs_chrom <- 
    HERVs_DE_anotados %>%
    filter(HERVs %in% c(modules$modules$`1`)) |>
    dplyr::select("Chrom","Start","End", "log2FoldChange", HERVs) %>%
    arrange(Chrom)

HERVs_chrom <- data.frame(
    chr = HERVs_chrom$Chrom,
    start = HERVs_chrom$Start,
    end = HERVs_chrom$End,
    value=HERVs_chrom$log2FoldChange,
    gene = HERVs_chrom$HERVs
    )

HERVs_chrom_up <-
  HERVs_chrom %>%
  filter(value>0)

color_control_alpha <- adjustcolor(color_control, alpha.f = 0.7)
color_covid_alpha   <- adjustcolor(color_covid, alpha.f = 1)
circos.initializeWithIdeogram(
  species = "hg38",
  ideogram.height = 0.1,
  track.height = 0.1,
  labels.cex = 1.5  
)

circos.genomicDensity(DEGs_chrom_list, col = c(color_covid_alpha, color_control_alpha), track.height = 0.3)
circos.genomicDensity(HERVs_chrom_up, col = colores[13], track.height = 0.2)

circos.clear()

### M2

DEG_chrom_FC <- 
  merge(mybdd[,c("chromosome_name","start_position","end_position", "hgnc_symbol")], genes_DE_anotados, by.x="hgnc_symbol", by.y="SYMBOL")

DEG_chrom <- DEG_chrom_FC %>%
  filter(!grepl("^H", chromosome_name)) %>%
  mutate(chr=paste("chr",chromosome_name, sep="")) %>%
  dplyr::select("chr","start_position","end_position", "log2FoldChange","hgnc_symbol", "ENTREZID") %>%
  arrange(chr)


DEG_chrom_sort <- 
  data.frame(
  chr = DEG_chrom$chr,
  start = DEG_chrom$start_position,
  end = DEG_chrom$end_position,
  value = DEG_chrom$log2FoldChange,
  gene = DEG_chrom$hgnc_symbol,
  gene_entrez = DEG_chrom$ENTREZID
)


DEG_chrom_sort_up <-
  DEG_chrom_sort %>%
    filter(gene_entrez %in% c(modules$modules$`2`)) |>
    dplyr::select(-gene_entrez) |>
  filter(value>0)

DEG_chrom_sort_down <-
  DEG_chrom_sort %>%
    filter(gene_entrez %in% c(modules$modules$`2`)) |>
    dplyr::select(-gene_entrez) |>
  filter(value<0)

DEGs_chrom_list = list(DEG_chrom_sort_up, DEG_chrom_sort_down)

HERVs_chrom <- 
    HERVs_DE_anotados %>%
    filter(HERVs %in% c(modules$modules$`2`)) |>
    dplyr::select("Chrom","Start","End", "log2FoldChange", HERVs) %>%
    arrange(Chrom)

HERVs_chrom <- data.frame(
    chr = HERVs_chrom$Chrom,
    start = HERVs_chrom$Start,
    end = HERVs_chrom$End,
    value=HERVs_chrom$log2FoldChange,
    gene = HERVs_chrom$HERVs
    )

HERVs_chrom_up <-
  HERVs_chrom %>%
  filter(value>0)

color_control_alpha <- adjustcolor(color_control, alpha.f = 0.7)
color_covid_alpha   <- adjustcolor(color_covid, alpha.f = 1)
circos.initializeWithIdeogram(
  species = "hg38",
  ideogram.height = 0.1,
  track.height = 0.1,
  labels.cex = 1.5 
)

circos.genomicDensity(DEGs_chrom_list, col = c(color_covid_alpha, color_control_alpha), track.height = 0.3)
circos.genomicDensity(HERVs_chrom_up, col = colores[13], track.height = 0.2)

circos.clear()

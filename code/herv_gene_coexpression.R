dds_prefil_all <-readRDS("data/dds_HERVs_filtered.rds")
dds_DE_all <- DESeq(dds_prefil_all)
res_all <- results(dds_DE_all) 

ddsSalmon_filtered_GC <- readRDS("data/dds_DEG_filtered_GC.rds")
ddsSalmon_filtered_GC <- DESeq(ddsSalmon_filtered_GC)
res_DEG_uni_GC <- results(ddsSalmon_filtered_GC)

HERVs_DE_ann_all <- readRDS("RDSs/HERVs_DE_ann_all.rds")
HERVs_DE_anotados <- HERVs_DE_ann_all |>
    filter(Class == "HERV")

up_hervs <- HERVs_DE_anotados |> 
    filter(log2FoldChange > 0)

#normalizar matriz conteos HERVs 
HERVs_mcount_vst <- vst(assay(dds_DE_all)[rownames(assay(dds_DE_all)) %in% up_hervs$HERVs, ])
#normalizar matriz conteos Genes
Genes_mcount_vst <- vst(assay(ddsSalmon_filtered_GC))

Genes_HERVs_mcounts_vst <- 
  bind_rows(
  as.data.frame(HERVs_mcount_vst), 
  as.data.frame(Genes_mcount_vst))
  
#comprobamos que se han unido bien las bbdd
Genes_HERVs_mcounts_vst["ERV316A3_10p13b", "COV1C"] == HERVs_mcount_vst["ERV316A3_10p13b", "COV1C"]
Genes_HERVs_mcounts_vst["1", "COV1C"] == Genes_mcount_vst["1", "COV1C"]

#formatting: genes as columns and samples as rows.
Genes_HERVs_mcounts_vst_t <- as.data.frame(t(Genes_HERVs_mcounts_vst))

Genes_HERVs_mcounts_vst_t_filtered <- filter_RNA_seq((filter_low_var(Genes_HERVs_mcounts_vst_t, pct = 0.7, type = "median")) , min_count=9, method = "mean")

threads_to_use <- 4
net_DE <- build_net(Genes_HERVs_mcounts_vst_t_filtered,
                    n_threads = threads_to_use, 
                    fit_cut_off=0.8, 
                    cor_func = "spearman",
                    keep_matrices = "cor"
                    )

modules <- detect_modules(Genes_HERVs_mcounts_vst_t_filtered, 
                            net_DE$network, 
                            detailled_result = TRUE,
                            merge_threshold = 0.25)

patients_data_all_GWENA <- readRDS("patients_data.rds)
patients_data_all_GWENA$condition <- as.character(patients_data_all_GWENA$condition)


patients_data_all_GWENA_sex <- patients_data_all_cov_RIN[rownames(Genes_HERVs_mcounts_vst_t_filtered),c("expediente", "condition", "sexo")]
patients_data_all_GWENA_sex$condition <- as.character(patients_data_all_GWENA_sex$condition)
patients_data_all_GWENA_sex$sexo <- as.character(patients_data_all_GWENA_sex$sexo)

phenotype_association <- associate_phenotype(
  modules$modules_eigengenes, 
  (patients_data_all_GWENA_sex %>% dplyr::select(c(condition, sexo))))

ME <- modules$modules_eigengenes

kMEs <- data.frame()

for (i in 1:length(modules$modules)){
    
    expr_module <- as.matrix(Genes_HERVs_mcounts_vst_t_filtered[, modules$modules[[i]]])
    
    common_samples <- intersect(rownames(expr_module), rownames(ME))
    expr_module <- expr_module[common_samples, , drop = FALSE]
    
    ME_mod <- ME[common_samples, i]
    kME <- cor(
        expr_module,
        ME_mod,
        method = "spearman",
        use = "pairwise.complete.obs"
        )
    
    kME_pvalue <- apply(
      expr_module,
      2,
      function(x) {
        cor.test(
          x,
          ME_mod,
          method = "spearman",
          exact = FALSE
        )$p.value
      }
    )
    
    kME_df <- data.frame(
        node = rownames(kME),
        kME = as.numeric(kME[, 1]),
        pvalue = kME_pvalue,
        padj = p.adjust(kME_pvalue, method = "BH"),
        module = i-1,
        stringsAsFactors = FALSE
        ) |> 
        filter(node %in%  modules$modules[[i]] & node %in% HERVs_DE_anotados$HERVs)
    
    if (nrow(kMEs) == 0) {
        kMEs <- kME_df
    } else {
        kMEs <- rbind(kMEs, kME_df)
    }
    
}


enriquecimientos <- list()

lista_enriquecer <- modules$modules$`3`[modules$modules$`3` %in% genes_DE_anotados$ENTREZID]
    
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

module_example <- modules$modules$`3`
g <- build_graph_from_sq_mat(net_DE$network[module_example, module_example])

herv_pattern <- HERVs_DE_anotados$HERVs

# =========================
# 3. GO terms de interés
# =========================

go_terms <-   list(
  inflammation = "GO:0006954",
  
  inflammasome = c(
    "GO:0061702",  # canonical inflammasome complex
    "GO:0097169",  # AIM2 inflammasome
    "GO:0140632",  # inflammasome assembly
    "GO:0044546",  # NLRP3 inflammasome assembly
    "GO:1900225",  # regulation NLRP3
    "GO:0070269",  # pyroptosis
    "GO:0050729"   # positive regulation inflammatory response
  ),
  
  chromatin = c(
    "GO:0016568",  # chromatin modification
    "GO:0006338",  # chromatin remodeling
    "GO:0006306"   # DNA methylation
  )
)

get_go_entrez <- function(go_ids) {
  
  genes <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = go_ids,
    keytype = "GOALL",
    columns = c("ENTREZID")
  ) %>%
    dplyr::pull(ENTREZID) %>%
    unique() %>%
    na.omit() %>%
    as.character()
  
  return(genes)
}

go_genes <- lapply(go_terms, get_go_entrez)

# =========================
# 4. Conversión Entrez -> símbolo
# =========================

all_go_entrez <- unique(unlist(go_genes))

entrez_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_go_entrez,
  keytype = "ENTREZID",
  columns = c("SYMBOL")
) %>%
  distinct(ENTREZID, .keep_all = TRUE) %>%
  mutate(
    ENTREZID = as.character(ENTREZID),
    SYMBOL = ifelse(is.na(SYMBOL), ENTREZID, SYMBOL)
  )

entrez_to_symbol <- setNames(entrez_map$SYMBOL, entrez_map$ENTREZID)

# =========================
# 5. Tabla de nodos y aristas
# =========================

vertex_df <- data.frame(
  node = V(g)$name,
  degree = degree(g),
  strength = strength(g, weights = E(g)$weight),
  stringsAsFactors = FALSE
) %>%
  mutate(
    type = ifelse(node %in% herv_pattern, "HERV", "GENE")
  )

edge_df <- as_data_frame(g, what = "edges") %>%
  mutate(abs_weight = abs(weight))

# =========================
# 6. Función para construir subgrafo por término
# =========================

build_term_graph <- function(term_name, edge_df, go_genes, herv_pattern, entrez_to_symbol,
                             top_hervs = 15, top_genes_per_herv = 10, lista_hervs = NULL) {
  
  term_gene_ids <- go_genes[[term_name]]
  
  # conservar solo aristas HERV-GENE del término
  term_edges <- edge_df %>%
    mutate(
      from_type = ifelse(from %in% herv_pattern, "HERV", "GENE"),
      to_type   = ifelse(to %in% herv_pattern, "HERV", "GENE")
    ) %>%
    filter(
      (from_type == "HERV" & to_type == "GENE" & to %in% term_gene_ids) |
      (from_type == "GENE" & to_type == "HERV" & from %in% term_gene_ids)
    ) %>%
    mutate(
      HERV = ifelse(from_type == "HERV", from, to),
      GENE = ifelse(from_type == "GENE", from, to)
    )
  
  if (nrow(term_edges) == 0) {
    message(paste("No hay conexiones para", term_name))
    return(NULL)
  }
  
  # priorizar HERVs por fuerza total dentro de ese término
  top_herv_names <- term_edges %>%
    group_by(HERV) %>%
    summarise(term_strength = sum(abs_weight, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(term_strength)) %>%
    slice_head(n = top_hervs) %>%
    pull(HERV)
  
    if (!is_null(lista_hervs)) {
        top_herv_names <- top_herv_names[top_herv_names %in% lista_hervs]
    }
  
  # quedarnos con esos HERVs
  term_edges <- term_edges %>%
    filter(HERV %in% top_herv_names)
  
  # por cada HERV, conservar solo sus top genes del término
  term_edges_small <- term_edges %>%
    group_by(HERV) %>%
    arrange(desc(abs_weight), .by_group = TRUE) %>%
    slice_head(n = top_genes_per_herv) %>%
    ungroup()
  
  # preparar tabla de nodos
  node_names <- unique(c(term_edges_small$HERV, term_edges_small$GENE))
  
  nodes_df <- data.frame(name = node_names, stringsAsFactors = FALSE) %>%
    mutate(
      type = ifelse(name %in% herv_pattern, "HERV", "GENE"),
      label = ifelse(type == "GENE",
                     ifelse(name %in% names(entrez_to_symbol), entrez_to_symbol[name], name),
                     name)
    )
  
  # construir grafo
  g_term <- graph_from_data_frame(
    term_edges_small %>% dplyr::select(from = HERV, to = GENE, weight),
    directed = FALSE,
    vertices = nodes_df
  )
  
  # atributos visuales
  V(g_term)$type_custom <- nodes_df$type[match(V(g_term)$name, nodes_df$name)]
  V(g_term)$label_custom <- nodes_df$label[match(V(g_term)$name, nodes_df$name)]
  V(g_term)$deg <- degree(g_term)
  V(g_term)$strg <- strength(g_term, weights = E(g_term)$weight)
  
  V(g_term)$color <- ifelse(V(g_term)$type_custom == "HERV", "tomato", "skyblue")
  V(g_term)$shape <- ifelse(V(g_term)$type_custom == "HERV", "square", "circle")
  V(g_term)$size <- ifelse(
    V(g_term)$type_custom == "HERV",
    12 + log1p(V(g_term)$strg + 1) * 3,
    6 + log1p(V(g_term)$deg + 1) * 2
  )
  
  E(g_term)$width <- 1 + 4 * abs(E(g_term)$weight) / max(abs(E(g_term)$weight), na.rm = TRUE)
  
  return(g_term)
}

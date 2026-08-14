data_samplequality_novogene <- read_csv("../data/HERVs/data_samplequality_novogene.csv")

quality_RNA <- data_samplequality_novogene |>
    dplyr::select(c( "Sample Name", "Concentration(ng/ul)", "Integrity value", "Sample QC Results")) |>
    mutate(group=sub("COV[0-9]*", "", data_samplequality_novogene$`Sample Name`))

for (i in 1:length(names(quality_RNA))) {
    if (names(quality_RNA)[i]=="Sample Name"){
        names(quality_RNA)[i]<-"expediente"}
    if (names(quality_RNA)[i]=="Concentration(ng/ul)"){
        names(quality_RNA)[i]<-"concentration_ngul"}
    if (names(quality_RNA)[i]=="Integrity value"){
        names(quality_RNA)[i]<-"RIN"}
    if (names(quality_RNA)[i]=="Sample QC Results"){
        names(quality_RNA)[i]<-"QC_cat"}
}

RNAseq_raw_quality <- read.delim("../data/HERVs/RNAseq_raw_quality.txt", header=TRUE)

RNAseq_raw_quality$Group <- NA
for (i in 1:dim(RNAseq_raw_quality)[1]) {
    if (sum(grep("COV[0-9]+C",RNAseq_raw_quality$Sample[i]))) 
        RNAseq_raw_quality$Group[i] <- "COVID"
    else if (sum(grep("COV[0-9]+E",RNAseq_raw_quality$Sample[i]))) 
        RNAseq_raw_quality$Group[i] <- "Control"
  
}

rownames(RNAseq_raw_quality) <- RNAseq_raw_quality$Sample

qualityanalysis <- data.frame(status=factor(), statistics=character(), file=character())

dirs_fastaqc <- list.dirs(path = "../data/HERVs/quality_analysis", recursive = FALSE,full.names = FALSE)

for (i in dirs_fastaqc) {
    newdb <- read.delim(paste("../data/HERVs/quality_analysis/",i,"/summary.txt",sep = ""), 
                      sep="\t",
                      col.names=c("status", "statistics", "file"))
    qualityanalysis <- rbind(newdb, qualityanalysis)
}

qualityanalysis$status <- factor(qualityanalysis$status)
qualityanalysis$file <- factor(qualityanalysis$file)

Bowtie_quality_raw <- read.table("../data/HERVs/multiqc_general_stats.txt",
                                 sep = "\t", 
                                 header = TRUE, 
                                 quote = "")

Run_stats_Telescope <- read.delim("../data/HERVs/Telescope_may2024/Run_stats_Telescope_encabezado2.txt")

files <- list.files(path = "../data/HERVs/Telescope_may2024/", 
                    pattern = "-TE_counts.tsv", 
                    all.files = FALSE,full.names = FALSE)

file_names <- sub("-TE_counts.tsv","", files)
matrixcount_all <- data.frame(transcript=character())


for (i in seq(length(file_names))) {
    file_bam <- read.delim(paste("../data/HERVs/Telescope_may2024/",files[i], sep ="")) 
    names(file_bam)[2] <-  file_names[i] 
    matrixcount_all <- merge(matrixcount_all,file_bam,by="transcript",all=TRUE) 
}


rownames(matrixcount_all) <- matrixcount_all$transcript
colnames(matrixcount_all) <- sub("_retro","",colnames(matrixcount_all))


RNA_seq_HERVs <- as.data.frame(matrixcount_all[,-1])
RNA_seq_HERVs["HERVs_NA",] <- colSums(is.na(RNA_seq_HERVs))
RNA_seq_HERVs["HERVs_0",] <- colSums(RNA_seq_HERVs==0, na.rm=TRUE)
RNA_seq_HERVs["HERVs_notO",] <- colSums(RNA_seq_HERVs>0, na.rm=TRUE)
RNA_seq_HERVs <- as.data.frame(t(RNA_seq_HERVs[c("__no_feature","HERVs_0","HERVs_notO","HERVs_NA"),]))

RNA_seq_HERVs_group <- merge(RNAseq_raw_quality[,c("Sample","Group"  )],RNA_seq_HERVs,by='row.names')

matrixcount_all <- matrixcount_all[rownames(matrixcount_all)!="__no_feature",]

patients_data_all <- read.table("../data/HERVs/groups_patients_all.txt")
patients_data_all$condition <- factor(patients_data_all$condition, levels=c("Control", "COVID"))

matrixcount_complete <- matrixcount_all[,rownames(patients_data_all)]

matrixcount_complete[is.na(matrixcount_complete)] <- 0

matrixcount_complete <- matrixcount_complete[,colnames(matrixcount_complete) != "COV20E"]



patient_HERvs_list <- as.data.frame(read_excel("../data/HERVs/patient_HERvs_list.xlsx"))

rownames(patient_HERvs_list) <- patient_HERvs_list$expediente

patients_data_all_cov <- merge(patients_data_all, patient_HERvs_list[,c("sexo", "edad", "expediente")], by='row.names', all=T)

rownames(patients_data_all_cov) <- patients_data_all_cov$expediente
patients_data_all_cov$sexo <- factor(patients_data_all_cov$sexo)

patients_data_all_cov_RIN <- merge(patients_data_all_cov, quality_RNA[,c("expediente","RIN", "QC_cat" )], by="expediente")
rownames(patients_data_all_cov_RIN) <- patients_data_all_cov_RIN$expediente
patients_data_all_cov_RIN$QC_cat <- factor(patients_data_all_cov_RIN$QC_cat)

patients_data_all_cov_RIN <- patients_data_all_cov_RIN[colnames(matrixcount_complete),]

patients_data_all_cov_RIN <- patients_data_all_cov_RIN |> filter(expediente != "COV20E")

patients_data_all_cov_RIN$RIN_centered <- patients_data_all_cov_RIN$RIN - mean(
    patients_data_all_cov_RIN$RIN,
    na.rm = TRUE
  )


temp <- read.delim("../data/Salmon_results_gcBias/COV1C_quant_gc/quant.sf")

Tx2Gene <- AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg38.knownGene, keys = as.vector(temp[, 1]),
    keytype = "TXNAME", columns = c("GENEID", "TXNAME"))
Tx2Gene <- Tx2Gene[!is.na(Tx2Gene$GENEID), ]


salmonQ_GC <- dir("../data/Salmon_results_gcBias/", 
               recursive = T, 
               pattern = "quant.sf", 
               full.names = T)

names(salmonQ_GC) <- sub("_quant_gc/quant.sf","", sub("../data/Salmon_results_gcBias//","", salmonQ_GC))



salmonCounts_GC <- tximport(salmonQ_GC, #We must provide the paths to Salmon
                         type = "salmon", #type of files
                         tx2gene = Tx2Gene) #data.frame of transcript to gene mapping

salmonCounts_GC$abundance <- salmonCounts_GC$abundance[,colnames(salmonCounts_GC$abundance) != "COV20E"]
salmonCounts_GC$counts <- salmonCounts_GC$counts[,colnames(salmonCounts_GC$counts) != "COV20E"]
salmonCounts_GC$length <- salmonCounts_GC$length[,colnames(salmonCounts_GC$length) != "COV20E"]


colnames(salmonCounts_GC$abundance) <- gsub("/", "", colnames(salmonCounts_GC$abundance))
colnames(salmonCounts_GC$counts) <- gsub("/", "", colnames(salmonCounts_GC$counts))
patients_data_Salmon_GC <- patients_data_all_cov_RIN[colnames(salmonCounts_GC$abundance),]

ddsSalmon_import_GC  <- DESeqDataSetFromTximport(salmonCounts_GC, 
                                      colData = patients_data_Salmon_GC, 
                                      design = ~ RIN_centered + condition)



keep_v02_S_GC <- rowSums(counts(ddsSalmon_import_GC)>10) >= 18
ddsSalmon_filtered_GC <- ddsSalmon_import_GC[keep_v02_S_GC,]

saveRDS(ddsSalmon_filtered_GC,"../data/dds_DEG_filtered_GC.rds")

patients_data_all_cov_RIN <- patients_data_all_cov_RIN |> 
    mutate(
        Exitus = c("No", "No","No","No","No","No","No","No","No","No",
                   "Yes", "Yes", "Yes", "Yes", "Yes", "No","No","No")
    ) |> 
    mutate(Exitus = factor(Exitus, levels = c("No", "Yes"), labels = c("Biopsy samples", "Autopsy samples")))


Bowtie_quality <- Bowtie_quality_raw |> 
    filter(!is.na(QualiMap_mqc_generalstats_qualimap_total_reads)) |>
    mutate(Sample_s=gsub('_sort_qm',"",Sample)) |> 
    dplyr::select(Sample_s, QualiMap_mqc_generalstats_qualimap_total_reads,QualiMap_mqc_generalstats_qualimap_mapped_reads,QualiMap_mqc_generalstats_qualimap_percentage_aligned, QualiMap_mqc_generalstats_qualimap_general_error_rate) |> 
    rename(Sample=Sample_s)


qc_table <- patients_data_all_cov_RIN |>
    dplyr::select(expediente, condition, sexo, edad, RIN, QC_cat, Exitus) |>
    left_join(
        Bowtie_quality |>
            dplyr::select(
                Sample,
                QualiMap_mqc_generalstats_qualimap_total_reads,
                QualiMap_mqc_generalstats_qualimap_mapped_reads,
                QualiMap_mqc_generalstats_qualimap_percentage_aligned,
                QualiMap_mqc_generalstats_qualimap_general_error_rate
                ),
        by = c("expediente" = "Sample")
    )

gene_detected <- sapply(salmonQ_GC, function(f) {
    x <- read.delim(f)
    sum(x$NumReads > 0, na.rm = TRUE)
    })


gene_detected_df <- data.frame(
    expediente = names(gene_detected),
    genes_detected = as.numeric(gene_detected)
    )

qc_table <- qc_table |>
    left_join(gene_detected_df, by = "expediente")

vars <- c(
    "RIN",
    "QualiMap_mqc_generalstats_qualimap_total_reads",
    "QualiMap_mqc_generalstats_qualimap_mapped_reads",
    "QualiMap_mqc_generalstats_qualimap_percentage_aligned",
    "QualiMap_mqc_generalstats_qualimap_general_error_rate",
    "genes_detected"
    )


# BIN-LEVEL ANALYSIS
#
# This script produces as output mutation score and amplification frequency for
# each chromosome (for each tumor type)

suppressMessages({
  library(readxl)
  library(tibble)
  library(ggplot2)
  library(tidyr)
  library(readr)
  library(stringr)
  library(data.table)
  library(dplyr)
  library(ggExtra)
  library(crayon)
  library(parallel)
})

cat("\n\n > This script calculates mu score and amplification frequency for each 1Mbp bin \n\n")
cat(
  " The computation will be performed for 23 cancer types and for \n several gene and mutation types (it will take a while)"
)

setwd("../")

tumor_types <- c(
  "BRCA",
  "LUAD",
  "LUSC",
  "CESC",
  "THCA",
  "HNSC",
  "PAAD",
  "COADREAD",
  "GBMLGG",
  # # #
  "SKCM",
  "BLCA",
  "PCPG",
  "PRAD",
  "KIRC",
  "MESO",
  "TGCT",
  "KIRP",
  "SARC",
  "LIHC",
  "ESCA",
  "STAD",
  "UCS",
  "OV"
)
    
mutation_types <- c(
  "all_mutations"
  
  ,"aggregation_causing"
  ,"non_aggregation_causing"
  
  ,"expressed_no0"
  ,"non_expressed_no0"

  ,"synonymous"
  ,"non_synonymous"
  # "silent"
  
  # "missense"
  ,"remove_OG"
  ,"remove_TSG"
  ,"remove_BOTH"
  
  ,"polyphen_highlyDamaging"
  ,"polyphen_moderatelyDamaging"
  ,"CADD_highlyDamaging_phred"
  ,"CADD_moderatelyDamaging_phred"
  
  ,"haploinsufficient"
  ,"non_haploinsufficient"
  ,"haploinsufficient_GHIS"
  ,"non_haploinsufficient_GHIS"
  
  # "haploinsufficient_damaging"
  # "haploinsufficient_nondamaging"
  # "non_haploinsufficient_damaging"
  
  # "haploinsufficient_synonymous"
  # "non_haploinsufficient_synonymous"
  # "haploinsufficient_non_synonymous"
  # "non_haploinsufficient_non_synonymous"
)

cat("\n\n Gene and mutation types are:\n\n")
print(mutation_types)

fixed_bin_length <- 1000000 # segmentation length (set at 1 Mbp)

# 5 10 20 30
segment_cutoff <- 20
cat("\n\n segment cutoff:", segment_cutoff, "\n\n")
stringent_mutations <- F

cores <- 23
# cores <- length(tumor_types)-8

# produce HAPLOINSUFFICIENT_GHIS score table
# ghis <- as.data.frame(readxl::read_xlsx("data/misc/nar-03716-met-n-2014-File006.xlsx", sheet = 3))
# colnames(ghis) <- c("ENSEMBL", "GHIS")
# library(org.Hs.eg.db)
# genes <- as.vector(ghis[,1])
# annots <- select(org.Hs.eg.db, keys=genes,
#                  columns="SYMBOL", keytype="ENSEMBL")
# result <- merge(ghis, annots, by.x="ENSEMBL", by.y="ENSEMBL")
# colnames(result) <- c("ENSEMBL", "GHIS", "Hugo_Symbol")
# write.table(result, file = "data/misc/ghis_scores.tsv", sep = "\t", quote = F, row.names = F)


for (mutation_type in mutation_types) {
  # compute mutations even in amplified and deleted regions
  if(mutation_type == "all_mutations_mutsWithinCNA"){
    mut_withinCNA <- T
  }else{
    mut_withinCNA <- F
  }
  
  # set/create the result folder according to mutation_type
  results_table_path <-
    paste0("results/tables/01_binLevelAnalysis/",
           mutation_type,
           "_vs_CN",ifelse(segment_cutoff=="20","",paste0("_0",segment_cutoff)),"/")                                                             ###### ADDED 025 HERE!!!!
  system(paste0("mkdir -p ", results_table_path))
  
  ## >> 1st Loop: cancer types << ----
  # for (tumor_type in tumor_types) {
  mclapply(tumor_types, mc.cores = cores, function(tumor_type){
    ## 1st Loop: load sCNAs and SNVs files ----
    cat(" \n\n\n >>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<< \n")
    cat(" >     ----------", paste0(tumor_type, " ---------- \n"))
    
    # Load SNVs file
    cat("\n > Load", tumor_type, "Mutation File \n")
    snv <-
      read.csv(file = paste0("data/FireBrowse_SNVs/", tumor_type, "_mutations.csv"))
    snv <-
      snv[!duplicated(snv[, c(6, 7, 8, 9, 17)]),] # remove duplicated mutations in the same patient
    # cat("  done!")
    
    # Load CNAs file
    cat("\n > Load", tumor_type, "SCNA File \n")
    scna <-
      read.delim(
        paste0(
          "data/FireBrowse_CNAs/",
          tumor_type,
          "_hg19_FireBrowse_totalCopyNumber_classifiedALL.tsv"
        ),
        header = T
      )
    # keep patients with both snv and scna data
    snv_filt <-
      snv[snv$patient_id %in% levels(factor(scna$patient_id)),]
    # keep patients with both snv and scna data
    scna <-
      scna[scna$patient_id %in% snv$patient_id,]
    
    common_patients <-
      levels(as.factor(snv[snv$patient_id %in% scna$patient_id,]$patient_id))
    # keep only non copy-neutral segments (amplified or deleted)
    scna <-
      scna[scna$Segment_Mean >= as.numeric(segment_cutoff)/100 | 
             scna$Segment_Mean <= -as.numeric(segment_cutoff)/100,]
    
    # cat("\n  done!")
    # cat("\n > common", tumor_type, "patients:", length(common_patients))
    
    ## 1st Loop: filter for gene properties ----
    # cat("\n > Filter for gene properties \n")
    if (mutation_type == "all_mutations") {
      cat(" no gene filtering")
    }
    
    if (mutation_type == "expressed_no0") {
      cat("\n Filter for gene type: expressed \n")
      tpm <- readRDS(file = "data/TCGA_tpm/PANCANCER_tpm_mean.rds.gz")
      tpm <- tpm[[tumor_type]]
      tpm <- tpm[tpm$median != 0,]
      snv_filt <-
        snv_filt[snv_filt$Hugo_Symbol %in% tpm$`Approved symbol`,]
      cat("Removed", dim(snv)[1] - dim(snv_filt)[1], "mutations")
    }
    if (mutation_type == "non_expressed_no0") {
      cat("\n Filter for gene type: expressed \n")
      tpm <- readRDS(file = "data/TCGA_tpm/PANCANCER_tpm_mean.rds.gz")
      tpm <- tpm[[tumor_type]]
      tpm <- tpm[tpm$median == 0,]
      snv_filt <-
        snv_filt[snv_filt$Hugo_Symbol %in% tpm$`Approved symbol`,]
      cat("Removed", dim(snv)[1] - dim(snv_filt)[1], "mutations")
    }
    
    if (mutation_type == "remove_OG") {
      cat("\n Filter for gene type: OGs \n")
      OGs <-
        read.csv(file = "data/CancerGenes/OG_list.tsv", sep = "\t")
      snv_filt <-
        snv_filt[!snv_filt$SWISSPROT %in%  OGs$Entry.name,]
      cat("Removed", dim(snv)[1] - dim(snv_filt)[1], "mutations")
      rm(OGs)
    }
    if (mutation_type == "remove_TSG") {
      cat("\n Filter for gene type: TSGs \n")
      TSGs <-
        read.csv(file = "data/CancerGenes/TSG_list.tab", sep = "\t")
      snv_filt <-
        snv_filt[!snv_filt$SWISSPROT %in%  TSGs$Entry.name,]
      cat("Removed", dim(snv)[1] - dim(snv_filt)[1], "mutations")
      rm(TSGs)
    }
    if (mutation_type == "remove_BOTH") {
      cat("\n Filter for gene type: both TSGs and OGs \n")
      OGs <-
        read.csv(file = "data/CancerGenes/OG_list.tsv", sep = "\t")
      TSGs <-
        read.csv(file = "data/CancerGenes/TSG_list.tab", sep = "\t")
      snv_filt <-
        snv_filt[!snv_filt$SWISSPROT %in%  TSGs$Entry.name,]
      snv_filt <-
        snv_filt[!snv_filt$SWISSPROT %in%  OGs$Entry.name,]
      cat("Removed", dim(snv)[1] - dim(snv_filt)[1], "mutations")
      rm(OGs)
      rm(TSGs)
    }
    
    
    if (mutation_type == "haploinsufficient_GHIS") {
      cat("\n Filter for haploinsufficient GHIS genes")
      ghis_scores <- read.table("data/misc/ghis_scores.tsv")
      colnames(ghis_scores) <- ghis_scores[1,]
      ghis_scores <- ghis_scores[-1,]
      ghis_scores$GHIS <- as.numeric(ghis_scores$GHIS)
      
      ghis_scores_intolerant <- ghis_scores[ghis_scores$GHIS >= 0.5,]
      snv_filt <-
        snv_filt[snv_filt$Hugo_Symbol %in% ghis_scores_intolerant$Hugo_Symbol,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
      snv_filt <- unique(snv_filt)
    }
    if (mutation_type == "non_haploinsufficient_GHIS") {
      cat("\n Filter for haploinsufficient GHIS genes")
      ghis_scores <- read.table("data/misc/ghis_scores.tsv")
      colnames(ghis_scores) <- ghis_scores[1,]
      ghis_scores <- ghis_scores[-1,]
      ghis_scores$GHIS <- as.numeric(ghis_scores$GHIS)
      
      ghis_scores_tolerant <- ghis_scores[ghis_scores$GHIS < 0.5,]
      snv_filt <-
        snv_filt[snv_filt$Hugo_Symbol %in% ghis_scores_tolerant$Hugo_Symbol,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
      snv_filt <- unique(snv_filt)
    }
    
    
    if (mutation_type == "haploinsufficient" | 
        mutation_type == "haploinsufficient_synonymous" |
        mutation_type == "haploinsufficient_non_synonymous") {
      cat("\n Filter for haploinsufficient genes")
      pLI_scores <-
        readxl::read_xlsx("data/misc/pLI_scores.xlsx", sheet = 2)
      pLI_scores <-
        pLI_scores[!is.na(pLI_scores$chr),] # remove genes within X and Y chromosomes
      pLI_scores_intolerant <- pLI_scores[pLI_scores$pLI >= 0.2,]
      snv_filt <-
        snv_filt[snv_filt$Hugo_Symbol %in% pLI_scores_intolerant$gene,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
      snv_filt <- unique(snv_filt)
    }
    if (mutation_type == "non_haploinsufficient"| 
        mutation_type == "non_haploinsufficient_synonymous" |
        mutation_type == "non_haploinsufficient__non_synonymous") {
      cat("\n Filter for non-haploinsufficient genes")
      pLI_scores <-
        readxl::read_xlsx("data/misc/pLI_scores.xlsx", sheet = 2)
      pLI_scores <-
        pLI_scores[!is.na(pLI_scores$chr),] # remove genes within X and Y chromosomes
      pLI_scores_tolerant <- pLI_scores[pLI_scores$pLI < 0.2,]
      snv_filt <-
        snv_filt[snv_filt$Hugo_Symbol %in% pLI_scores_tolerant$gene,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
      snv_filt <- unique(snv_filt)
    }
    
    if (mutation_type == "haploinsufficient_damaging") {
      # cat("\n Filter for haploinsufficient_damaging genes")
      pLI_scores <-
        readxl::read_xlsx("data/misc/pLI_scores.xlsx", sheet = 2)
      pLI_scores <-
        pLI_scores[!is.na(pLI_scores$chr),] # remove genes within X and Y chromosomes
      pLI_scores_intolerant <- pLI_scores[pLI_scores$pLI >= 0.2,]
      snv_filt <-
        snv_filt[snv_filt$Hugo_Symbol %in% pLI_scores_intolerant$gene,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
      snv_filt <- unique(snv_filt)
      # damaging CADD
      snv_filt <- snv_filt[as.numeric(snv_filt$CADD_phred) >= 20,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
    }
    if (mutation_type == "haploinsufficient_nondamaging") {
      # cat("\n Filter for haploinsufficient_nondamaging genes")
      pLI_scores <-
        readxl::read_xlsx("data/misc/pLI_scores.xlsx", sheet = 2)
      pLI_scores <-
        pLI_scores[!is.na(pLI_scores$chr),] # remove genes within X and Y chromosomes
      pLI_scores_intolerant <- pLI_scores[pLI_scores$pLI >= 0.2,]
      snv_filt <-
        snv_filt[snv_filt$Hugo_Symbol %in% pLI_scores_intolerant$gene,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
      snv_filt <- unique(snv_filt)
      # non-damaging CADD 
      snv_filt <- snv_filt[as.numeric(snv_filt$CADD_raw) < 3.5,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
    }
    if (mutation_type == "non_haploinsufficient_damaging") {
      # cat("\n Filter for haploinsufficient_nondamaging genes")
      pLI_scores <-
        readxl::read_xlsx("data/misc/pLI_scores.xlsx", sheet = 2)
      pLI_scores <-
        pLI_scores[!is.na(pLI_scores$chr),] # remove genes within X and Y chromosomes
      pLI_scores_intolerant <- pLI_scores[pLI_scores$pLI < 0.2,]
      snv_filt <-
        snv_filt[snv_filt$Hugo_Symbol %in% pLI_scores_intolerant$gene,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
      snv_filt <- unique(snv_filt)
      # non-damaging CADD 
      snv_filt <- snv_filt[as.numeric(snv_filt$CADD_phred) >= 20,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
    }
    
    
    ## 1st Loop: filter for mutation properties ----
    # cat("\n > Filter for Mutations Type \n")
    if (mutation_type == "all_mutations") {
      cat("\n no mutation filtering\n")
    }
    
    if (mutation_type == "polyphen_highlyDamaging") {
      # cat("\n Filter for polyphen_highlyDamaging mutations")
      mean(as.numeric(snv_filt$Polyphen2_HVAR_score), na.rm = T)
      snv_filt <- snv_filt[!is.na(snv_filt$Polyphen2_HVAR_score),]
      snv_filt <-
        snv_filt[as.numeric(snv_filt$Polyphen2_HVAR_score) >= 0.6,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
      snv_filt <- unique(snv_filt)
    }
    if (mutation_type == "polyphen_moderatelyDamaging") {
      # cat("\n Filter for polyphen_moderatelyDamaging mutations")
      snv_filt <- snv_filt[!is.na(snv_filt$Polyphen2_HVAR_score),]
      snv_filt <-
        snv_filt[as.numeric(snv_filt$Polyphen2_HVAR_score) <= 0.3,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
      snv_filt <- unique(snv_filt)
    }
    
    if (mutation_type == "CADD_highlyDamaging") {
      # cat("\n Filter for CADD_highlyDamaging mutations")
      snv_filt <- snv_filt[as.numeric(snv_filt$CADD_phred) >= 20,]
      snv_filt <- snv_filt[as.numeric(snv_filt$CADD_raw) >= 3.5,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
    }
    if (mutation_type == "CADD_highlyDamaging_phred") {
      snv_filt <- snv_filt[as.numeric(snv_filt$CADD_phred) >= 20,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
    }
    if (mutation_type == "CADD_moderatelyDamaging") {
      # cat("\n Filter for CADD_moderatelyDamaging mutations")
      snv_filt <- snv_filt[as.numeric(snv_filt$CADD_phred) < 20,]
      snv_filt <- snv_filt[as.numeric(snv_filt$CADD_raw) < 3.5,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
    }
    if (mutation_type == "CADD_moderatelyDamaging_phred") {
      # cat("\n Filter for CADD_moderatelyDamaging mutations")
      snv_filt <- snv_filt[as.numeric(snv_filt$CADD_phred) < 20,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
    }
    
    if (mutation_type == "synonymous" | 
        mutation_type == "haploinsufficient_synonymous" |
        mutation_type == "non_haploinsufficient_synonymous") {
      print("filter for syn")
      snv_filt <- snv_filt[snv_filt$Consequence == "synonymous_variant",]
      # snv_filt <- snv_filt[snv_filt$Variant_Classification == "Silent",]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
    }
    if (mutation_type == "non_synonymous" | 
        mutation_type == "haploinsufficient_non_synonymous" |
        mutation_type == "non_haploinsufficient_non_synonymous") {
      print("filter for non_syn")
      snv_filt <- snv_filt[snv_filt$Consequence != "synonymous_variant",]
      # snv_filt <- snv_filt[snv_filt$Consequence == "missense_variant",]
      # snv_filt <- snv_filt[snv_filt$Variant_Classification == "Missense_Mutation",]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
    }
    
    
    if (mutation_type == "aggregation_causing") {
      # cat("\n Filter for aggregation_causing mutations")
      snv_filt <-
        snv_filt[snv_filt$Aggregation > 5000 &
                   snv_filt$FoldChange > 1,]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
      snv_filt <- unique(snv_filt)
      print(paste0("Aggregating mutations: ", dim(snv_filt)))
    }
    if (mutation_type == "non_aggregation_causing") {
      # cat("\n Filter for non_aggregation_causing mutations")
      snv_filt <-
        snv_filt[!(snv_filt$Aggregation > 5000 &
                     snv_filt$FoldChange > 1),]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
      snv_filt <- unique(snv_filt)
      print(paste0("Non-aggregating mutations: ", dim(snv_filt)))
    }
    
    # filter for mutations Variant_Classification
    if (mutation_type == "nonsense") {
      snv_filt <-
        snv_filt[snv_filt$Variant_Classification == "Nonsense_Mutation",]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
    }
    if (mutation_type == "silent") {
      snv_filt <- snv_filt[snv_filt$Variant_Classification == "Silent",]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
    }
    if (mutation_type == "missense") {
      cat(" Filter for non_aggregation_causing mutations")
      
      snv_filt <- snv_filt[snv_filt$Consequence == "missense_variant" | 
                             snv_filt$Variant_Classification == "Missense_Mutation",]
      snv_filt <- snv_filt[!is.na(snv_filt$Chromosome),]
    }
    
    if (mutation_type == "randomized_10000") {
      snv_filt <- sample_n(snv, 10000)
    }
    
    rm(snv)
    gc(full = T)
    
    
    ##  initialize table for chromosome and arm amplification frequencies ----
    cat("\n\n >> Start", tumor_type, "chromosome Loop \n")
    freq_ampl_chr <- t(data.frame("p", "p+q", "q"))
    colnames(freq_ampl_chr) <- "arms"
    rownames(freq_ampl_chr) <- NULL
    freq_ampl_chr <- as.data.frame(freq_ampl_chr)
    
    ## >> 2nd Loop: chromosome << ----
    for (chr in 1:22) {
      gc(full = TRUE)
      # print(chr)
      
      ## 2nd loop: load chromosome region lengths ----
      chr_info <-
        read.table("data/misc/chr_info_h19.txt", header = TRUE)
      chr_arms <-
        read.table(file = "data/misc/cytoBand.txt", header = T)
      chr_arms[, 2:3] <-
        apply(chr_arms[, 2:3] / fixed_bin_length, 2, as.integer)
      
      all_patients <-
        levels(as.factor(as.character(scna$patient_id)))
      
      ## 2nd loop: create temporary tables for SNVs and CNAs ----
      temp_cna <- scna[scna$Chromosome == chr,]
      temp_snv <- snv_filt[snv_filt$Chromosome == chr,]
      
      ## 2nd loop: define binning parameters ----
      n_bins <-
        as.integer(chr_info[chr_info$Chromosome == paste0("chr", chr),]$Length /
                     fixed_bin_length)
      length_bin <-
        as.integer(chr_info[chr_info$Chromosome == paste0("chr", chr),]$Length /
                     n_bins)
      
      ## 2nd loop: compute start and length for each CNA according to bins ----
      end <- length_bin
      start_bin <- NULL
      temp_cna_length <- NULL
      i <- 1
      
      # this function assign each CNA event to segments 
      whichbin <- function(data, end, i) {
        if (any(data$Start < end) == TRUE) {
          temp2 <- data[data$Start < end,]
          temp2 <-
            cbind(
              temp2,
              start_bin = rep(i, nrow(temp2)),
              bin_length = (temp2$End - temp2$Start) / fixed_bin_length
            )
        }
        return(temp2)
      }
      
      # apply the whichbin() function
      for (i in 1:n_bins) {
        if (nrow(temp_cna) == 0) {
          break
        }
        if (any(temp_cna$Start < end)) {
          temp_cna_length <-
            rbind(temp_cna_length, whichbin(temp_cna, end, i))
        }
        temp_cna <- temp_cna[!temp_cna$Start < end,]
        end <- end + length_bin
      }
      
      # keep CNA events higher that 0,5 Mbp (half of bin size)
      temp_cna_length <-
        temp_cna_length[temp_cna_length$length > fixed_bin_length / 2,]
      amplified_patients <-
        levels(factor(temp_cna_length$patient_id))
      
      if(stringent_mutations != T){
        # produce the table
        write.table(
          amplified_patients,
          file = paste0(
            results_table_path,
            tumor_type,
            "_amplified_pts_chr",
            chr,
            ".txt"
          )
        )
      }
      
      ## 2nd loop: compute chromosome and arm amplification frequencies ----
      temp_cna_length_backup <- temp_cna_length
      temp_cna_length <-
        temp_cna_length %>% filter(Segment_Mean >= as.numeric(segment_cutoff)/100)
      
      if(dim(temp_cna_length)[1] == 0){
        ampl <-
          cbind(rbind("p", "p+q", "q", "chr"),
                rbind(
                  p = NA,
                  q = NA,
                  'p+q' = NA,
                  chr = chr
                ))
        colnames(ampl) <- c("arms", "Freq")
        ampl <- as.data.frame(ampl)
        freq_ampl_chr <-
          full_join(freq_ampl_chr, ampl,  by = "arms")
      }else{
        arm <- data.frame()
        for (i in 1:nrow(temp_cna_length)) {
          p <-
            chr_arms[chr_arms$chromosme == paste0("chr", chr),]$start[1]:chr_arms[chr_arms$chromosme == paste0("chr", chr),]$end[1] #p
          q <-
            chr_arms[chr_arms$chromosme == paste0("chr", chr),]$start[2]:chr_arms[chr_arms$chromosme == paste0("chr", chr),]$end[2] #q
          
          if (temp_cna_length[i,]$classified != "chromosomal") {
            cna <-
              temp_cna_length[i,]$start_bin:(temp_cna_length[i,]$start_bin + temp_cna_length[i,]$bin_length)
            
            if (sum(cna %in% p) > sum(cna %in% q)) {
              arm <-
                rbind(
                  arm,
                  cbind(
                    patient_id = temp_cna_length[i,]$patient_id,
                    classified = temp_cna_length[i,]$classified,
                    arm = "p"
                  )
                )
            } else{
              arm <-
                rbind(
                  arm,
                  cbind(
                    patient_id = temp_cna_length[i,]$patient_id,
                    classified = temp_cna_length[i,]$classified,
                    arm = "q"
                  )
                )
            }
          } else{
            arm <-
              rbind(
                arm,
                cbind(
                  patient_id = temp_cna_length[i,]$patient_id,
                  classified = temp_cna_length[i,]$classified,
                  arm = "p+q"
                )
              )
          }
        }
        
        temp_cna_length <- cbind(temp_cna_length, arm = arm$arm)
        
        temp_cna_length[temp_cna_length$classified == "arm" &
                          temp_cna_length$arm == "q" &
                          temp_cna_length$start_bin + temp_cna_length$bin_length > chr_arms[chr_arms$chromosme == paste0("chr", chr),]$end[1] &
                          temp_cna_length$start_bin < 0.5 * chr_arms[chr_arms$chromosme == paste0("chr", chr),]$end[1],]$arm <-
          rep("p+q", length(temp_cna_length[temp_cna_length$classified == "arm" &
                                              temp_cna_length$arm == "q" &
                                              temp_cna_length$start_bin +
                                              temp_cna_length$bin_length > chr_arms[chr_arms$chromosme == paste0("chr", chr),]$end[1] &
                                              temp_cna_length$start_bin < 0.5 *
                                              chr_arms[chr_arms$chromosme == paste0("chr", chr),]$end[1],]$arm))
        
        ampl <-
          as.data.frame(rbind(as.data.frame(
            table(temp_cna_length[temp_cna_length$Segment_Mean > 0 &
                                    temp_cna_length$classified != "focal",]$arm) / length(common_patients)
          )))
        
        if (nrow(ampl) != 0) {
          colnames(ampl) <- c("arms", "Freq")
          ampl <-
            rbind(c(arms = c("chr", chr), Freq = as.numeric(chr)), ampl)
          freq_ampl_chr <-
            full_join(freq_ampl_chr, ampl,  by = "arms")
        } else{
          ampl <-
            cbind(rbind("p", "p+q", "q", "chr"),
                  rbind(
                    p = NA,
                    q = NA,
                    'p+q' = NA,
                    chr = chr
                  ))
          colnames(ampl) <- c("arms", "Freq")
          ampl <- as.data.frame(ampl)
          freq_ampl_chr <-
            full_join(freq_ampl_chr, ampl,  by = "arms")
        }
      }
      
      
      ## 2nd loop: load segmented chromosome/gene structure ----
      bin_gene <-
        read.table(
          paste0(
            "data/ChromosomeGeneStructure/chr_",
            chr,
            "_binSize_",
            fixed_bin_length,
            ".txt"
          )
        )
      
      ## 2nd loop: compute amplification and deletion frequencies ----
      chr_bins <- data.frame(bins = 1:n_bins)
      
      temp_cna_length <- temp_cna_length_backup
      
      for (i in 1:nrow(temp_cna_length)) {
        if (temp_cna_length$bin_length[i] > 0.5) {
          y <-
            temp_cna_length$bin_length[i] - as.integer(temp_cna_length$bin_length[i])
          if (y >= 0.5) {
            x <-
              data.frame(
                bins = temp_cna_length$start_bin[i]:as.integer(
                  temp_cna_length$start_bin[i] + temp_cna_length$bin_length[i]
                ),
                c(rep(
                  1, as.integer(temp_cna_length$bin_length[i] + 1)
                ))
              )
          } else{
            x <-
              data.frame(
                bins = temp_cna_length$start_bin[i]:as.integer(
                  temp_cna_length$start_bin[i] + temp_cna_length$bin_length[i] - 1
                ),
                c(rep(
                  1, as.integer(temp_cna_length$bin_length[i])
                ))
              )
          }
        } else{
          next
        }
        
        chr_bins <- full_join(chr_bins, x, by = "bins")
        
        colnames(chr_bins) <-
          c(colnames(chr_bins)[-length(colnames(chr_bins))], paste0("CNA", i))
      }
      rm(i)
      rm(x)
      rm(y)
      
      chr_bins <- t(as.matrix(chr_bins))
      chr_bins <- as.data.frame(chr_bins)
      chr_bins <- chr_bins[-1,]
      
      # CNA shorter than half of the bin size (0.5 Mbp in this case) were removed from calculation
      chr_bins <-
        cbind(
          Segment_Mean = temp_cna_length[temp_cna_length$bin_length > 0.5,]$Segment_Mean,
          classified = temp_cna_length[temp_cna_length$bin_length > 0.5,]$classified,
          patient = temp_cna_length[temp_cna_length$bin_length > 0.5,]$patient_id,
          cna_length_del = temp_cna_length[temp_cna_length$bin_length > 0.5 |
                                             temp_cna_length$Segment_Mean <= -as.numeric(segment_cutoff)/100,]$length,
          cna_length_ampl = temp_cna_length[temp_cna_length$bin_length > 0.5 |
                                              temp_cna_length$Segment_Mean >= as.numeric(segment_cutoff)/100,]$length,
          chr_bins
        )
      
      cna_freq_ampl <-
        as.numeric(colSums(chr_bins[chr_bins$Segment_Mean >= as.numeric(segment_cutoff)/100, 6:ncol(chr_bins)], na.rm = TRUE) /
                     length(common_patients))
      cna_freq_del <-
        as.numeric(colSums(chr_bins[chr_bins$Segment_Mean <= -as.numeric(segment_cutoff)/100, 6:ncol(chr_bins)], na.rm = TRUE) /
                     length(common_patients))
      cna_freq_total <-
        as.numeric(colSums(chr_bins[, 6:ncol(chr_bins)], na.rm = TRUE) / length(common_patients))
      
      bin_gene <- cbind(
        bin_gene,
        cna_freq_ampl = cna_freq_ampl,
        cna_freq_del = cna_freq_del,
        cna_freq_total = cna_freq_total
      )
      
      if(mut_withinCNA == T){
        unamplified_patients <- all_patients
      }else{
        unamplified_patients <-
          all_patients[!(all_patients %in% chr_bins$patient)]
      }
      
      ## 2nd loop: compute the mutation score in copy-neutral regions ----
      # (normalized for the number of patients (1), the coding region (2) and thus log10 (3))
      chr_bins[is.na(chr_bins)] <- 0
      
      chr_bins_pt <- data.frame()
      for (pt in levels(factor(chr_bins$patient))) {
        colSums(chr_bins[chr_bins$patient == pt,][,-c(1:5)])
        chr_bins_pt <-
          rbind(chr_bins_pt, c(pt, as.numeric(colSums(chr_bins[chr_bins$patient == pt,][,-c(1:5)]))))
      }
      colnames(chr_bins_pt) <- c("patients", 1:n_bins)
      
      x <- NULL
      end <- length_bin
      start_bin <- 0
      mutations_raw <- NULL
      mutations_norm <- NULL
      n_patients <- NULL
      
      for (i in 1:n_bins) {
        if(mut_withinCNA == T){
          n_pts <- length(all_patients)
          mut <- temp_snv[str_sub(temp_snv$Tumor_Sample_Barcode, end = 12) %in% all_patients,]
          mut <-
            mut[mut$Start_Position >= start_bin &
                  mut$End_Position < end,]
          mut_a <- as.data.frame(table(mut$patient_id))        
        }else{
          n_pts <-
            length(c(chr_bins_pt[chr_bins_pt[, i + 1] == 0,]$patient, unamplified_patients))
          mut <-
            temp_snv[str_sub(temp_snv$Tumor_Sample_Barcode, end = 12) %in% 
                       c(chr_bins_pt[chr_bins_pt[, i +1] == 0,]$patient, 
                         unamplified_patients),]
          mut <-
            mut[mut$Start_Position >= start_bin &
                  mut$End_Position < end,]
          mut_a <- as.data.frame(table(mut$patient_id))        
        }
        
        if (dim(mut_a)[1] != 0) {
          colnames(mut_a) <- c("patient_id", "mutations_raw")
          x <- sum(mut_a$mutations_raw , na.rm = T)
          mutations_raw <- c(mutations_raw, x)
          mutations_norm <- c(mutations_norm, x / n_pts)
        } else{
          mutations_raw <- c(mutations_raw, 0)
          mutations_norm <- c(mutations_norm, 0)
        }
        n_patients <- c(n_patients, n_pts)
        start_bin <- start_bin + length_bin
        end <- end + length_bin
      }
      
      bin_gene_mut <-
        cbind(bin_gene, mutations_raw, mutations_norm, n_patients)
      bin_gene_mut <- as.data.frame(bin_gene_mut)
      
      ## end of the 2nd loop: write chromosome table ----
      write.table(
        bin_gene_mut,
        file = paste0(
          results_table_path,
          tumor_type,
          "_chr",
          chr,
          "_",
          as.integer(fixed_bin_length),
          "BIN_table.txt"
        )
      )
      
    }
    
    ## end of the 1st loop: write chromosome/arm amplification frequencies ----
    colnames(freq_ampl_chr) <- NULL
    freq_ampl_chr[is.na(freq_ampl_chr)] <- 0
    freq_ampl_chr <- rbind(freq_ampl_chr, chr = as.numeric(freq_ampl_chr[4,])+as.numeric(freq_ampl_chr[5,]))
    freq_ampl_chr <- freq_ampl_chr[-c(4,5),]
    if(stringent_mutations != T){
      write.table(freq_ampl_chr,
                  file = paste0(results_table_path, tumor_type, 
                                "_chrAmpFreq",
                                ifelse(segment_cutoff=="20","",paste0("_0",segment_cutoff))
                                ,".txt"))                               
    }  
    cat("\n\n> >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Analysis on", tumor_type, " ended.\n\n")
  })
  # }
}

cat(
  "\n\n OUTPUT of the script: \n \t (1) raw tables path: results/tables/01_binLevelAnalysis/ \n"
)

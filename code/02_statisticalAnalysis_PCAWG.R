# SEGMENT-LEVEL ANALYSIS
#
# This script produces as output mutation score and amplification frequency for
# each chromosome (for each tumor type) with different bin sizes and conditions

suppressMessages({
  require(ggplot2)
  require(stringr)
  library(readr)
  require(mgcv)
  require(ggpubr)
  library(optparse)
  library(dplyr)
  library(utils)
})

# set optparse parameters 
option_list = list(
  make_option(
    c("-t", "--tables"),
    type = "character",
    default = "y",
    help = "Options are: ([y]/n)
                It requires some minutes (up to 1 h). 
                It produces bin level correlation tables for different segment
                lengths and conditions.",
    metavar = ""
  ),
  make_option(
    c("-s", "--statistics"),
    type = "character",
    default = "y",
    help = "Options are: ([y]/n)
                It produces statistics and plots.",
    metavar = ""
  )
)

opt_parser = OptionParser(option_list = option_list)

opt = parse_args(opt_parser)

if (opt$tables == "y" & opt$statistics == "y") {
  cat("\n\n >> You chose default options: \n\t(1) --tables 'y';\n\t(2) --statistics 'y' \n\n")
}

if (!any(opt$tables %in% c("y", "n"))) {
  print_help(opt_parser)
  stop("typo in the analysis flag, plase see above for the available options!",
       call. = FALSE)
}

cat(
  "\n\n > This script \n\n\t (1) produce correlations between amplification frequency and mu score for different bin sizes:
    \t\t - from 1 to 50 Mbp
    \t\t - chromosome-arm level
    \t\t - entire chromosome level \n
    \n\t (2) produces correlations and figures \n\n\n"
)

produce_tables <- TRUE
if (opt$tables == "n") {
  produce_tables <- FALSE
}
produce_statistics <- TRUE
if (opt$statistics == "n") {
  produce_statistics <- FALSE
}

setwd("../")

tumor_types <- paste0(c(
  "Breast",
  "Brain",
  "Kidney",
  "Liver",
  "Prostate",
  "Ovary",
  "Skin",
  "Pancreas"
), ".merged")

# tumor_types <- levels(factor(columns$V1))
# tumor_types <- as.data.frame(cbind(do.call(rbind, str_split(tumor_types, "-")), tumor_types))
# tumor_types <- tumor_types$tumor_types
# tumor_types <- tumor_types[!tumor_types %in% c("CLLE-ES", "CMDI-UK", "LAML-KR", "THCA-US")]

intron_intergenic <- T

if (F) {
  columns <- read.table("data/PCAWG/patient_per_tumortype.tsv",
                        skip = 1)
  tumor_types <- levels(factor(columns$V1))
  tumor_types <- as.data.frame(cbind(do.call(rbind, str_split(tumor_types, "-")), tumor_types))
  tumor_types <- tumor_types$tumor_types
  tumor_types <- tumor_types[!(tumor_types %in% c("BLCA-US", "BRCA-US", "CESC-US",
                                                  "COAD-US", "DLBC-US", "GBM-US",
                                                  "HNSC-US", "KICH-US", "KIRC-US",
                                                  "KIRP-US", "LAML-KR", "LGG-US", 
                                                  "LIHC-US", "LUAD-US", "LUSC-US",
                                                  "OV-US", "PRAD-US", "READ-US",
                                                  "SARC-US", "SKCM-US", "STAD-US",
                                                  "THCA-US", "UCEC-US"
  ))]
}

n_mutations <- data.frame()
if (produce_tables) {
  # output table directory
  results_table_path <-
    paste0("results/tables/02_produceStatistics_PCAWG/")
  system(paste0("mkdir -p ", results_table_path))
  
  # initialize parameters
  parameters_chr <- data.frame()
  parameters_arm <- data.frame()
  parameters <- data.frame()
  tierfinal <- data.frame()
  
  # segmentation lengths (default is 36Mbp)
  segment_lengths <- c(36)
  
  ## >> 1st loop (optional): varying segment lengths << ----
  for (segment_length in segment_lengths) {
    # 36 was a choose as default (best) cutoff (Fig.1d)
    cat("\n > bin size (Mbp): ", segment_length, "\n")
    
    if (segment_length != 36) {
      conditions <- c(
        "amplifications"
        ,"deletions")
    } else{
      conditions <- c(
        "amplifications"
        ,"deletions"
        ,"coding"
        ,"noncoding"
      )
    }
    
    ## >> 2nd loop: mutation/gene type conditions << ----
    source_table_path <-
      paste0("results/tables/01_binLevel_PCAWG/")
    
    ## >> 3rd loop: CANCER TYPE << ----
    suppressMessages({
      for (tumor_type in tumor_types) {
        # 3rd loop: load chromosome tables (output from script 01_mainAnalysis.R) ----
        chr1 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr1_1000000BIN_table.txt"
            )
          ), chr = 1)
        chr2 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr2_1000000BIN_table.txt"
            )
          ), chr = 2)
        chr3 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr3_1000000BIN_table.txt"
            )
          ), chr = 3)
        chr4 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr4_1000000BIN_table.txt"
            )
          ), chr = 4)
        chr5 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr5_1000000BIN_table.txt"
            )
          ), chr = 5)
        chr6 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr6_1000000BIN_table.txt"
            )
          ), chr = 6)
        chr7 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr7_1000000BIN_table.txt"
            )
          ), chr = 7)
        chr8 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr8_1000000BIN_table.txt"
            )
          ), chr = 8)
        chr9 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr9_1000000BIN_table.txt"
            )
          ), chr = 9)
        chr10 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr10_1000000BIN_table.txt"
            )
          ), chr = 10)
        chr11 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr11_1000000BIN_table.txt"
            )
          ), chr = 11)
        chr12 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr12_1000000BIN_table.txt"
            )
          ), chr = 12)
        chr13 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr13_1000000BIN_table.txt"
            )
          ), chr = 13)
        chr14 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr14_1000000BIN_table.txt"
            )
          ), chr = 14)
        chr15 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr15_1000000BIN_table.txt"
            )
          ), chr = 15)
        chr16 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr16_1000000BIN_table.txt"
            )
          ), chr = 16)
        chr17 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr17_1000000BIN_table.txt"
            )
          ), chr = 17)
        chr18 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr18_1000000BIN_table.txt"
            )
          ), chr = 18)
        chr19 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr19_1000000BIN_table.txt"
            )
          ), chr = 19)
        chr20 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr20_1000000BIN_table.txt"
            )
          ), chr = 20)
        chr21 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr21_1000000BIN_table.txt"
            )
          ), chr = 21)
        chr22 <-
          cbind(read.table(
            file = paste0(
              source_table_path,
              tumor_type,
              "_chr22_1000000BIN_table.txt"
            )
          ), chr = 22)
        
        
        # 3rd loop: add a new column ----
        # to create segments of desired segment_length and remove bins that do not contains genes
        add.col <- function(df) {
          # df <- df[df$gene_count != 0, ]
          df <- cbind(df,
                      resize = rep(
                        1:round(dim(df)[1] / segment_length + 0.5),
                        each = round(dim(df)[1] / round(
                          dim(df)[1] / segment_length + 0.5
                        ) + 0.5)
                      )[1:nrow(df)])
          return(df)
        }
        
        chr1 <-  add.col(chr1)
        chr2 <-  add.col(chr2)
        chr3 <-  add.col(chr3)
        chr4 <-  add.col(chr4)
        chr5 <-  add.col(chr5)
        chr6 <-  add.col(chr6)
        chr7 <-  add.col(chr7)
        chr8 <-  add.col(chr8)
        chr9 <-  add.col(chr9)
        chr10 <- add.col(chr10)
        chr11 <- add.col(chr11)
        chr12 <- add.col(chr12)
        chr13 <- add.col(chr13)
        chr14 <- add.col(chr14)
        chr15 <- add.col(chr15)
        chr16 <- add.col(chr16)
        chr17 <- add.col(chr17)
        chr18 <- add.col(chr18)
        chr19 <- add.col(chr19)
        chr20 <- add.col(chr20)
        chr21 <- add.col(chr21)
        chr22 <- add.col(chr22)
        
        
        # 3rd loop: create segments of desired segment_length ----
        merge.bins_coding <- function(df) {
          df2 <- data.frame()
          chr <- as.numeric(df$chr[1])
          for (i in levels(factor(df$resize))) {
            df <- df[!(df$gene_count <= 1 & df$length_noncoding <= 50000), ]
            
            start <- df[df$resize == i, ]$bin_start[1]
            end <-
              df[df$resize == i, ]$bin_end[length(df[df$resize == i, ]$bin_start)]
            gene_count <- sum(df[df$resize == i, ]$gene_count)
            
            length_coding <- sum(df[df$resize == i, ]$length_coding)
            length_noncoding <- sum(df[df$resize == i, ]$length_noncoding)
            
            cna_freq_ampl <- mean(df[df$resize == i, ]$cna_freq_ampl)
            cna_freq_del <-  mean(df[df$resize == i, ]$cna_freq_del)
            
            mutations_raw <-
              sum(as.numeric(df[df$resize == i, ]$mutations_raw))
            mutations_coding <-
              sum(as.numeric(df[df$resize == i, ]$mutations_coding_wintron))
            mutations_coding_wointron <-
              sum(as.numeric(df[df$resize == i, ]$mutations_coding_wointron))
            mutations_noncoding_intron <-
              sum(as.numeric(df[df$resize == i, ]$mutations_noncoding_intron))
            mutations_noncoding_intergenic <-
              sum(as.numeric(df[df$resize == i, ]$mutations_noncoding_intergenic))
            
            gene_id <-
              paste0(df[df$resize == i, ]$gene_id, collapse = "")
            
            if (dim(df[df$resize == i, ])[1] == 0) {
              next
            }
            if (mutations_raw == 0) {
              mutations_raw <- 0
              mutations_coding <- 0
              mutations_coding_wointron <- 0
              mutations_noncoding_intron <- 0
              mutations_noncoding_intergenic <- 0
              mutations_normPT <- 0
              mutations_coding_normPT <- 0
              mutations_coding_wointron_normPT <- 0
              mutations_noncoding_intron_normPT <- 0
              mutations_noncoding_intergenic_normPT <- 0
            }else{
              mutations_normPT <-
                mean(
                  as.numeric(df[df$resize == i, ]$mutations_norm),
                  na.rm = T
                )
              mutations_coding_normPT <-
                mean(
                  as.numeric(df[df$resize == i, ]$mutations_coding_norm),
                  na.rm = T
                )
              mutations_coding_wointron_normPT <-
                mean(
                  as.numeric(df[df$resize == i, ]$mutations_coding_wointron_norm),
                  na.rm = T
                )
              
              if(intron_intergenic == F){
                mutations_noncoding_intron <- 0
                mutations_noncoding_intergenic <- 0
                mutations_noncoding_intron_normPT <- 0
                mutations_noncoding_intergenic_normPT <- 0
              }else{                  
                mutations_noncoding_intron_normPT <-
                  mean(
                    as.numeric(df[df$resize == i, ]$mutations_noncoding_intron_norm),
                    na.rm = T
                  )
                mutations_noncoding_intergenic_normPT <-
                  mean(
                    as.numeric(df[df$resize == i, ]$mutations_noncoding_intergenic_norm),
                    na.rm = T
                  )
              }
              
              # }
            }
            
            # mutations are normalized according to Eq. 2
            
            df2 <- rbind.data.frame(
              df2,
              c(
                gene_count,
                start,
                end,
                length_coding,
                length_noncoding,
                cna_freq_ampl,
                cna_freq_del,
                mutations_raw,
                mutations_coding,
                mutations_coding_wointron,
                mutations_noncoding_intron,
                mutations_noncoding_intergenic,
                mutations_normPT,
                mutations_coding_normPT,
                mutations_coding_wointron_normPT,
                mutations_noncoding_intron_normPT,
                mutations_noncoding_intergenic_normPT,
                gene_id,
                gene_count,
                chr,
                resize = i
              ),
              stringsAsFactors = FALSE
            )
          }
          
          colnames(df2) <- c(
            "ene_count",
            "start",
            "end",
            "length_coding",
            "length_noncoding",
            "cna_freq_ampl",
            "cna_freq_del",
            "mutations_raw",
            "mutations_coding",
            "mutations_coding_wointron",
            "mutations_noncoding_intron",
            "mutations_noncoding_intergenic",
            "mutations_normPT",
            "mutations_coding_normPT",
            "mutations_coding_wointron_normPT",
            "mutations_noncoding_intron_normPT",
            "mutations_noncoding_intergenic_normPT",
            "gene_id",
            "gene_count",
            "chr",
            "resize"
          )
          
          return(df2)
        }
        merge.bins_noncoding <- function(df) {
          df2 <- data.frame()
          chr <- as.numeric(df$chr[1])
          for (i in levels(factor(df$resize))) {
            df <- df[!(df$gene_count <= 1 & df$length_noncoding <= 50000), ]
            
            start <- df[df$resize == i, ]$bin_start[1]
            end <-
              df[df$resize == i, ]$bin_end[length(df[df$resize == i, ]$bin_start)]
            gene_count <- sum(df[df$resize == i, ]$gene_count)
            
            length_coding <- sum(df[df$resize == i, ]$length_coding)
            length_noncoding <- sum(df[df$resize == i, ]$length_noncoding)
            
            cna_freq_ampl <- mean(df[df$resize == i, ]$cna_freq_ampl)
            cna_freq_del <-  mean(df[df$resize == i, ]$cna_freq_del)
            
            mutations_raw <-
              sum(as.numeric(df[df$resize == i, ]$mutations_raw))
            mutations_coding <-
              sum(as.numeric(df[df$resize == i, ]$mutations_coding_wintron))
            mutations_coding_wointron <-
              sum(as.numeric(df[df$resize == i, ]$mutations_coding_wointron))
            mutations_noncoding_intron <-
              sum(as.numeric(df[df$resize == i, ]$mutations_noncoding_intron))
            mutations_noncoding_intergenic <-
              sum(as.numeric(df[df$resize == i, ]$mutations_noncoding_intergenic))
            
            gene_id <-
              paste0(df[df$resize == i, ]$gene_id, collapse = "")
            
            if (dim(df[df$resize == i, ])[1] == 0) {
              next
            }
            if (mutations_raw == 0) {
              mutations_raw <- 0
              mutations_coding <- 0
              mutations_coding_wointron <- 0
              mutations_noncoding_intron <- 0
              mutations_noncoding_intergenic <- 0
              mutations_normPT <- 0
              mutations_coding_normPT <- 0
              mutations_coding_wointron_normPT <- 0
              mutations_noncoding_intron_normPT <- 0
              mutations_noncoding_intergenic_normPT <- 0
            }else{
              mutations_normPT <-
                mean(
                  as.numeric(df[df$resize == i, ]$mutations_norm),
                  na.rm = T
                )
              mutations_coding_normPT <-
                mean(
                  as.numeric(df[df$resize == i, ]$mutations_coding_norm),
                  na.rm = T
                )
              mutations_coding_wointron_normPT <-
                mean(
                  as.numeric(df[df$resize == i, ]$mutations_coding_wointron_norm),
                  na.rm = T
                )
              
              if(intron_intergenic == F){
                mutations_noncoding_intron <- 0
                mutations_noncoding_intergenic <- 0
                mutations_noncoding_intron_normPT <- 0
                mutations_noncoding_intergenic_normPT <- 0
              }else{                  
                mutations_noncoding_intron_normPT <-
                  mean(
                    as.numeric(df[df$resize == i, ]$mutations_noncoding_intron_norm),
                    na.rm = T
                  )
                mutations_noncoding_intergenic_normPT <-
                  mean(
                    as.numeric(df[df$resize == i, ]$mutations_noncoding_intergenic_norm),
                    na.rm = T
                  )
              }
              
              # }
            }
            
            # mutations are normalized according to Eq. 2
            
            df2 <- rbind.data.frame(
              df2,
              c(
                gene_count,
                start,
                end,
                length_coding,
                length_noncoding,
                cna_freq_ampl,
                cna_freq_del,
                mutations_raw,
                mutations_coding,
                mutations_coding_wointron,
                mutations_noncoding_intron,
                mutations_noncoding_intergenic,
                mutations_normPT,
                mutations_coding_normPT,
                mutations_coding_wointron_normPT,
                mutations_noncoding_intron_normPT,
                mutations_noncoding_intergenic_normPT,
                gene_id,
                gene_count,
                chr,
                resize = i
              ),
              stringsAsFactors = FALSE
            )
          }
          
          colnames(df2) <- c(
            "ene_count",
            "start",
            "end",
            "length_coding",
            "length_noncoding",
            "cna_freq_ampl",
            "cna_freq_del",
            "mutations_raw",
            "mutations_coding",
            "mutations_coding_wointron",
            "mutations_noncoding_intron",
            "mutations_noncoding_intergenic",
            "mutations_normPT",
            "mutations_coding_normPT",
            "mutations_coding_wointron_normPT",
            "mutations_noncoding_intron_normPT",
            "mutations_noncoding_intergenic_normPT",
            "gene_id",
            "gene_count",
            "chr",
            "resize"
          )
          
          return(df2)
        }
        
        chr1_coding <-  merge.bins_coding(chr1)
        chr2_coding <-  merge.bins_coding(chr2)
        chr3_coding <-  merge.bins_coding(chr3)
        chr4_coding <-  merge.bins_coding(chr4)
        chr5_coding <-  merge.bins_coding(chr5)
        chr6_coding <-  merge.bins_coding(chr6)
        chr7_coding <-  merge.bins_coding(chr7)
        chr8_coding <-  merge.bins_coding(chr8)
        chr9_coding <-  merge.bins_coding(chr9)
        chr10_coding <- merge.bins_coding(chr10)
        chr11_coding <- merge.bins_coding(chr11)
        chr12_coding <- merge.bins_coding(chr12)
        chr13_coding <- merge.bins_coding(chr13)
        chr14_coding <- merge.bins_coding(chr14)
        chr15_coding <- merge.bins_coding(chr15)
        chr16_coding <- merge.bins_coding(chr16)
        chr17_coding <- merge.bins_coding(chr17)
        chr18_coding <- merge.bins_coding(chr18)
        chr19_coding <- merge.bins_coding(chr19)
        chr20_coding <- merge.bins_coding(chr20)
        chr21_coding <- merge.bins_coding(chr21)
        chr22_coding <- merge.bins_coding(chr22)
        
        tier_coding <-
          rbind(
            chr1_coding,
            chr2_coding,
            chr3_coding,
            chr4_coding,
            chr5_coding,
            chr6_coding,
            chr7_coding,
            chr8_coding,
            chr9_coding,
            chr10_coding,
            chr11_coding,
            chr12_coding,
            chr13_coding,
            chr14_coding,
            chr15_coding,
            chr16_coding,
            chr17_coding,
            chr18_coding,
            chr19_coding,
            chr20_coding,
            chr21_coding,
            chr22_coding
          )
        
        tier_coding[, c(1:18, 20:21)] <-
          apply(tier_coding[, c(1:18, 20:21)], 2, as.numeric)
        
        
        chr1_noncoding <-  merge.bins_noncoding(chr1)
        chr2_noncoding <-  merge.bins_noncoding(chr2)
        chr3_noncoding <-  merge.bins_noncoding(chr3)
        chr4_noncoding <-  merge.bins_noncoding(chr4)
        chr5_noncoding <-  merge.bins_noncoding(chr5)
        chr6_noncoding <-  merge.bins_noncoding(chr6)
        chr7_noncoding <-  merge.bins_noncoding(chr7)
        chr8_noncoding <-  merge.bins_noncoding(chr8)
        chr9_noncoding <-  merge.bins_noncoding(chr9)
        chr10_noncoding <- merge.bins_noncoding(chr10)
        chr11_noncoding <- merge.bins_noncoding(chr11)
        chr12_noncoding <- merge.bins_noncoding(chr12)
        chr13_noncoding <- merge.bins_noncoding(chr13)
        chr14_noncoding <- merge.bins_noncoding(chr14)
        chr15_noncoding <- merge.bins_noncoding(chr15)
        chr16_noncoding <- merge.bins_noncoding(chr16)
        chr17_noncoding <- merge.bins_noncoding(chr17)
        chr18_noncoding <- merge.bins_noncoding(chr18)
        chr19_noncoding <- merge.bins_noncoding(chr19)
        chr20_noncoding <- merge.bins_noncoding(chr20)
        chr21_noncoding <- merge.bins_noncoding(chr21)
        chr22_noncoding <- merge.bins_noncoding(chr22)
        
        tier_noncoding <-
          rbind(
            chr1_noncoding,
            chr2_noncoding,
            chr3_noncoding,
            chr4_noncoding,
            chr5_noncoding,
            chr6_noncoding,
            chr7_noncoding,
            chr8_noncoding,
            chr9_noncoding,
            chr10_noncoding,
            chr11_noncoding,
            chr12_noncoding,
            chr13_noncoding,
            chr14_noncoding,
            chr15_noncoding,
            chr16_noncoding,
            chr17_noncoding,
            chr18_noncoding,
            chr19_noncoding,
            chr20_noncoding,
            chr21_noncoding,
            chr22_noncoding
          )
        
        tier_noncoding[, c(1:18, 20:21)] <-
          apply(tier_noncoding[, c(1:18, 20:21)], 2, as.numeric)
        
        # output table of the analysis
        write_tsv(
          tier_coding,
          file = paste0(
            results_table_path,
            tumor_type,
            "_",
            segment_length,
            "Mbp_table_coding.tsv"
          )
        )
        write_tsv(
          tier_noncoding,
          file = paste0(
            results_table_path,
            tumor_type,
            "_",
            segment_length,
            "Mbp_table_noncoding.tsv"
          )
        )
        
        
        n_mutations <- rbind(n_mutations,
                             cbind(tumor_type = tumor_type,
                                   coding = sum(tier_coding$mutations_coding, na.rm = T),
                                   non_coding = sum(tier_noncoding$mutations_noncoding, na.rm = T)))
        
        
          corP_amplifications <- try(cor.test(log10(tier_coding$mutations_normPT), 
                                              tier_coding$cna_freq_ampl, method = "pearson"))
          corS_amplifications <- try(cor.test(log10(tier_coding$mutations_normPT), 
                                              tier_coding$cna_freq_ampl, method = "spearman"))
          
          corP_deletions <- try(cor.test(log10(tier_coding$mutations_normPT), 
                                         tier_coding$cna_freq_del, method = "pearson"))
          corS_deletions <- try(cor.test(log10(tier_coding$mutations_normPT), 
                                         tier_coding$cna_freq_del, method = "spearman"))
          
          corP_coding <- try(cor.test(log10(tier_coding$mutations_coding_normPT), 
                                      tier_coding$cna_freq_ampl, method = "pearson"))
          corS_coding <- try(cor.test(log10(tier_coding$mutations_coding_normPT), 
                                      tier_coding$cna_freq_ampl, method = "spearman"))
          
          corP_coding_deletions <- try(cor.test(log10(tier_coding$mutations_coding_normPT), 
                                                tier_coding$cna_freq_del, method = "pearson"))
          corS_coding_deletions <- try(cor.test(log10(tier_coding$mutations_coding_normPT), 
                                                tier_coding$cna_freq_del, method = "spearman"))
          
          corP_coding_wointron <- try(cor.test(log10(tier_noncoding$mutations_coding_wointron_normPT), 
                                               tier_noncoding$cna_freq_ampl, method = "pearson"))
          corS_coding_wointron <- try(cor.test(log10(tier_noncoding$mutations_coding_wointron_normPT), 
                                               tier_noncoding$cna_freq_ampl, method = "spearman"))
          
          corP_coding_wointron_deletions <- try(cor.test(log10(tier_noncoding$mutations_coding_wointron_normPT), 
                                                         tier_noncoding$cna_freq_del, method = "pearson"))
          corS_coding_wointron_deletions <- try(cor.test(log10(tier_noncoding$mutations_coding_wointron_normPT), 
                                                         tier_noncoding$cna_freq_del, method = "spearman"))
          
          if(intron_intergenic){
            corP_noncoding_intron <- try(cor.test(log10(tier_noncoding$mutations_noncoding_intron_normPT), 
                                                  tier_noncoding$cna_freq_ampl, method = "pearson"))
            corS_noncoding_intron <- try(cor.test(log10(tier_noncoding$mutations_noncoding_intron_normPT), 
                                                  tier_noncoding$cna_freq_ampl, method = "spearman"))
            
            corP_noncoding_intergenic <- try(cor.test(log10(tier_noncoding$mutations_noncoding_intergenic_normPT), 
                                                      tier_noncoding$cna_freq_ampl, method = "pearson"))
            corS_noncoding_intergenic <- try(cor.test(log10(tier_noncoding$mutations_noncoding_intergenic_normPT), 
                                                      tier_noncoding$cna_freq_ampl, method = "spearman"))
            
            corP_noncoding_intron_deletions <- try(cor.test(log10(tier_noncoding$mutations_noncoding_intron_normPT), 
                                                            tier_noncoding$cna_freq_del, method = "pearson"))
            corS_noncoding_intron_deletions <- try(cor.test(log10(tier_noncoding$mutations_noncoding_intron_normPT), 
                                                            tier_noncoding$cna_freq_del, method = "spearman"))
            
            corP_noncoding_intergenic_deletions <- try(cor.test(log10(tier_noncoding$mutations_noncoding_intergenic_normPT), 
                                                                tier_noncoding$cna_freq_del, method = "pearson"))
            corS_noncoding_intergenic_deletions <- try(cor.test(log10(tier_noncoding$mutations_noncoding_intergenic_normPT), 
                                                                tier_noncoding$cna_freq_del, method = "spearman"))
          }
        
        tryCatch({
          parameters <- rbind(parameters,  rbind(c(
            tumor_type = tumor_type,
            segment_length = segment_length,
            condition = "amplifications",
            corP = corP_amplifications$estimate,
            p.corP = corP_amplifications$p.value,
            corS = corS_amplifications$estimate,
            p.corS = corS_amplifications$p.value
          ),c(
            tumor_type = tumor_type,
            segment_length = segment_length,
            condition = "coding",
            corP = corP_coding$estimate,
            p.corP = corP_coding$p.value,
            corS = corS_coding$estimate,
            p.corS = corS_coding$p.value
          ),c(
            tumor_type = tumor_type,
            segment_length = segment_length,
            condition = "noncoding_intron",
            corP = corP_noncoding_intron$estimate,
            p.corP = corP_noncoding_intron$p.value,
            corS = corS_noncoding_intron$estimate,
            p.corS = corS_noncoding_intron$p.value
          ),c(
            tumor_type = tumor_type,
            segment_length = segment_length,
            condition = "noncoding_intergenic",
            corP = corP_noncoding_intergenic$estimate,
            p.corP = corP_noncoding_intergenic$p.value,
            corS = corS_noncoding_intergenic$estimate,
            p.corS = corS_noncoding_intergenic$p.value
          )
          ))
        })
      }
    })
  }
}

colnames(parameters) <-
  c("tumorType",
    "segment_length",
    "condition",
    "corP",
    "p.corP",
    "corS",
    "p.corS")
parameters[, 4:7] <- apply(parameters[, 4:7], 2, as.numeric)


aggregat <-
  parameters[parameters$condition == "coding" |
               parameters$condition == "noncoding_intron" |
               parameters$condition == "noncoding_intergenic",]

p1 <-
  ggplot(aggregat, aes(
    x = as.factor(tumorType),
    y = corS,
    fill = condition
  )) +
  geom_bar(color = "black",
           stat = "identity",
           position = position_dodge()) +
  theme_classic() + scale_fill_brewer(palette = "Blues")
p2 <-
  ggplot(aggregat, aes(
    x = as.factor(tumorType),
    y = log10(p.corS),
    fill = condition
  )) +
  geom_bar(color = "black",
           stat = "identity",
           position = position_dodge()) +
  theme_classic() + scale_fill_brewer(palette = "Reds")

ggarrange(p1, p2, nrow = 2)

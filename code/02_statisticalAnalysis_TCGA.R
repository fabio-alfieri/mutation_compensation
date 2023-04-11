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
produce_statistics <- F
if (opt$statistics == "n") {
  produce_statistics <- FALSE
}

setwd("../")

tumor_types <- c(
  "LUAD",
  "LUSC",
  "BRCA",
  "CESC",
  "THCA",
  "HNSC",
  "PAAD",
  "COADREAD",
  "GBMLGG",
  # #
  "SKCM",
  "OV",
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
  "UCS"
)

segment_cutoffs <- c(20)
all_conditions <- data.frame()

if (produce_tables) {
  if (F) {
    # compute extremes uncovered by CN sequencing method (AffymetrixSNP6.0)
    averageChrCoverage <- data.frame()
    for (tumor_type in tumor_types[1:8]) {
      scna <-
        read.delim(
          paste0(
            "data/FireBrowse_CNAs/",
            tumor_type,
            "_hg19_FireBrowse_totalCopyNumber_classifiedALL.tsv"
          ),
          header = T
        )
      for (chr in 1:22) {
        min <- min(scna[scna$Chromosome == chr, ]$Start)
        max <- max(scna[scna$Chromosome == chr, ]$End)
        averageChrCoverage <-
          rbind(
            averageChrCoverage,
            cbind(
              tumor_type = tumor_type,
              chr = chr,
              start = min,
              end = max
            )
          )
      }
    }
    write.table(file = "data/FireBrowse_CNAs/averageChrCoverage.txt", averageChrCoverage[, c(2:4)])
  }
  extremes <- "remove" # "remove" o ""
  # extremes were removed because CNAs and SNVs comes from different
  # sequencing methods, thus chromosome extremes only contains SNVs, but no CNAs.
  averageChrCoverage <-
    read.table("data/FireBrowse_CNAs/averageChrCoverage.txt")
  averageChrCoverage <- averageChrCoverage[c(1:22), ]
  
  # initialize tables
  # output table directory
  results_table_path <-
    paste0("results/tables/02_produceStatistics/")
  system(paste0("mkdir -p ", results_table_path))
  # plot table directory
  results_plot_path <-
    paste0("results/plots/02_tumorCorrelations/")
  system(paste0("mkdir -p ", results_plot_path))
  
  # initialize parameters
  parameters_chr <- data.frame()
  parameters_arm <- data.frame()
  parameters <- data.frame()
  tierfinal <- data.frame()
  
  # segmentation lengths (default is 36Mbp)
  segment_lengths <- c(15,18,20,22,25,28,30,33,36)
  
  for(segment_cutoff_muts in segment_cutoffs_muts){
    cat("\n > MUTATION segment_mean cutoff: ", segment_cutoff_muts)
    
    for(segment_cutoff in segment_cutoffs){
      cat("\n >> (1) analyis: produce correlations at different bin sizes \n")
      
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
            
            ,"expressed_no0"
            ,"non_expressed_no0"

            ,"silent"
            ,"synonymous"
            ,"non_synonymous"
            
            ,"aggregation_causing"
            ,"non_aggregation_causing"

            ,"missense"
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
            ,"haploinsufficient_synonymous"
          )
        }
        
        if (segment_length == 1) {
          # produce chromosome/arm level correlations
          chr_arm <- "yes" # "" or "yes"
        } else{
          chr_arm <- ""
        }
        
        ## >> 2nd loop: mutation/gene type conditions << ----
        for (condition in conditions) {
          if(condition == "amplifications" | condition == "deletions"){
            source_table_path <-
              paste0("results/tables/01_binLevelAnalysis/all_mutations_vs_CN",
                     ifelse(segment_cutoff=="20","",paste0("_0",segment_cutoff)),"/")
          }else{
            source_table_path <-
              paste0("results/tables/01_binLevelAnalysis/",
                     condition,
                     "_vs_CN",ifelse(segment_cutoff=="20","",paste0("_0",segment_cutoff)),"/")
          }
          
          cat(" > correlation condition: ", condition, "\n")
          
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
              
              
              # 3rd loop: remove the extremes of the chromosomes (uncovered by CNAs) ----
              if (extremes == "remove") {
                remove.extremes <- function(df) {
                  df <-
                    df[!df$bin_end <= as.numeric(averageChrCoverage[averageChrCoverage$chr == df$chr[1], ]$start), ]
                  df <-
                    df[!df$bin_start >= as.numeric(averageChrCoverage[averageChrCoverage$chr == df$chr[1], ]$end), ]
                  return(df)
                }
                
                chr1 <-  remove.extremes(chr1)
                chr2 <-  remove.extremes(chr2)
                chr3 <-  remove.extremes(chr3)
                chr4 <-  remove.extremes(chr4)
                chr5 <-  remove.extremes(chr5)
                chr6 <-  remove.extremes(chr6)
                chr7 <-  remove.extremes(chr7)
                chr8 <-  remove.extremes(chr8)
                chr9 <-  remove.extremes(chr9)
                chr10 <- remove.extremes(chr10)
                chr11 <- remove.extremes(chr11)
                chr12 <- remove.extremes(chr12)
                chr13 <- remove.extremes(chr13)
                chr14 <- remove.extremes(chr14)
                chr15 <- remove.extremes(chr15)
                chr16 <- remove.extremes(chr16)
                chr17 <- remove.extremes(chr17)
                chr18 <- remove.extremes(chr18)
                chr19 <- remove.extremes(chr19)
                chr20 <- remove.extremes(chr20)
                chr21 <- remove.extremes(chr21)
                chr22 <- remove.extremes(chr22)
              }
              
              # 3rd loop: add a new column ----
              # to create segments of desired segment_length and remove bins that do not contains genes
              add.col <- function(df) {
                df <- df[df$gene_count != 0, ]
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
              merge.bins <- function(df) {
                df2 <- data.frame()
                chr <- as.numeric(df$chr[1])
                for (i in levels(factor(df$resize))) {
                  df <- df[!(df$gene_count <= 1 & df$length_perc <= 0.05), ]
                  start <- df[df$resize == i, ]$bin_start[1]
                  end <-
                    df[df$resize == i, ]$bin_end[length(df[df$resize == i, ]$bin_start)]
                  gene_count <- sum(df[df$resize == i, ]$gene_count)
                  length_coding <-
                    sum(df[df$resize == i, ]$length_coding)
                  length_perc <- mean(df[df$resize == i, ]$length_perc)
                  cna_freq_ampl <-
                    mean(df[df$resize == i, ]$cna_freq_ampl)
                  cna_freq_del <- mean(df[df$resize == i, ]$cna_freq_del)
                  mutations_raw <-
                    sum(as.numeric(df[df$resize == i, ]$mutations_raw))
                  gene_id <-
                    paste0(df[df$resize == i, ]$gene_id, collapse = "")
                  if (mutations_raw == 0) {
                    mutations_raw <- 0
                  }
                  if (dim(df[df$resize == i, ])[1] == 0) {
                    next
                  }
                  mutations_norm <-
                    mean(
                      as.numeric(df[df$resize == i, ]$mutations_norm) / as.numeric(df[df$resize == i, ]$length_coding),
                      na.rm = T
                    )
                  # mutations are normalized according to Eq. 2
                  
                  df2 <- rbind.data.frame(
                    df2,
                    c(
                      gene_count,
                      start,
                      end,
                      length_coding,
                      length_perc,
                      cna_freq_ampl,
                      cna_freq_del,
                      mutations_raw,
                      mutations_norm,
                      gene_id,
                      chr,
                      resize = i
                    ),
                    stringsAsFactors = FALSE
                  )
                }
                
                colnames(df2) <- c(
                  "gene_count",
                  "start",
                  "end",
                  "length_coding",
                  "length_perc",
                  "cna_freq_ampl",
                  "cna_freq_del",
                  "mutations_raw",
                  "mutations_norm",
                  "gene_id",
                  "chr",
                  "resize"
                )
                
                return(df2)
              }
              
              chr1 <-  merge.bins(chr1)
              chr2 <-  merge.bins(chr2)
              chr3 <-  merge.bins(chr3)
              chr4 <-  merge.bins(chr4)
              chr5 <-  merge.bins(chr5)
              chr6 <-  merge.bins(chr6)
              chr7 <-  merge.bins(chr7)
              chr8 <-  merge.bins(chr8)
              chr9 <-  merge.bins(chr9)
              chr10 <- merge.bins(chr10)
              chr11 <- merge.bins(chr11)
              chr12 <- merge.bins(chr12)
              chr13 <- merge.bins(chr13)
              chr14 <- merge.bins(chr14)
              chr15 <- merge.bins(chr15)
              chr16 <- merge.bins(chr16)
              chr17 <- merge.bins(chr17)
              chr18 <- merge.bins(chr18)
              chr19 <- merge.bins(chr19)
              chr20 <- merge.bins(chr20)
              chr21 <- merge.bins(chr21)
              chr22 <- merge.bins(chr22)
              
              tier <-
                rbind(
                  chr1,
                  chr2,
                  chr3,
                  chr4,
                  chr5,
                  chr6,
                  chr7,
                  chr8,
                  chr9,
                  chr10,
                  chr11,
                  chr12,
                  chr13,
                  chr14,
                  chr15,
                  chr16,
                  chr17,
                  chr18,
                  chr19,
                  chr20,
                  chr21,
                  chr22
                )
              
              tier[, c(1:9, 11:12)] <-
                apply(tier[, c(1:9, 11:12)], 2, as.numeric)
              
              all_conditions <- rbind(all_conditions,
                                      cbind(tier,
                                            condition = condition,
                                            segment_cutoff_muts = segment_cutoff_muts,
                                            tumor_type = tumor_type))
              
              # output table of the analysis
              write_tsv(
                tier,
                file = paste0(
                  results_table_path,
                  condition,
                  "_",
                  tumor_type,
                  "_",
                  segment_length,
                  "Mbp_table",
                  ifelse(segment_cutoff=="020","",paste0("_",segment_cutoff)),
                  ".tsv"
                )
              )
              
              # assess correlation at chromosome/arm level
              if (chr_arm == "yes") {
                ampl_freq <-
                  read.table(file = paste0(source_table_path, tumor_type, "_chrAmpFreq",
                                           ifelse(segment_cutoff=="020","",paste0("_",segment_cutoff)),
                                           ".txt"))
                chr_ampl_mut <- data.frame()
                for (chr in 1:22) {
                  mut <-
                    mean(tier[tier$chr == chr, ]$mutations_norm / tier[tier$chr == chr, ]$length_coding)
                  if (chr == 13 |
                      chr == 14 | chr == 15 | chr == 21 | chr == 22) {
                    ampl <- ampl_freq[3, chr + 2]
                  } else{
                    ampl <- ampl_freq[2, chr + 2]
                  }
                  chr_ampl_mut <-
                    rbind(chr_ampl_mut,
                          cbind(chr, ampl, mut),
                          stringsAsFactors = FALSE)
                }
                (corP_chr <-
                    cor.test(
                      chr_ampl_mut$ampl,
                      log10(chr_ampl_mut$mut),
                      method = "pearson"
                    ))
                (corS_chr <-
                    cor.test(
                      chr_ampl_mut$ampl,
                      log10(chr_ampl_mut$mut),
                      method = "spearman"
                    ))
                
                parameters_chr <-
                  rbind(
                    parameters_chr,
                    c(
                      tumor_type = tumor_type,
                      type = "chr",
                      condition = condition,
                      corP = corP_chr$estimate,
                      p.valP = corP_chr$p.value,
                      corS = corS_chr$estimate,
                      p.valS = corS_chr$p.value
                    ),
                    stringsAsFactors = FALSE
                  )
                
                cytoband <-
                  read.table(file = "data/misc/cytoBand.txt", header = T)
                arm_ampl_mut <- data.frame()
                for (chr in 1:22) {
                  l <- cytoband[cytoband$chromosme == paste0("chr", chr), ][1, ]$end
                  p <- tier[tier$chr == chr & tier$end <= l, ]
                  q <- tier[tier$chr == chr & tier$end > l, ]
                  mut_p <- mean(p$mutations_norm / p$length_coding)
                  mut_q <- mean(q$mutations_norm / q$length_coding)
                  if (chr == 13 |
                      chr == 14 | chr == 15 | chr == 21 | chr == 22) {
                    ampl_q <- ampl_freq[3, chr + 2]
                    arm_ampl_mut <-
                      rbind(
                        arm_ampl_mut,
                        cbind(
                          chr = paste0(chr, "q"),
                          mut = mut_q,
                          ampl = ampl_q
                        ),
                        stringsAsFactors = FALSE
                      )
                  } else{
                    ampl_p <- ampl_freq[2, chr + 2]
                    ampl_q <- ampl_freq[3, chr + 2]
                    arm_ampl_mut <-
                      rbind(
                        arm_ampl_mut,
                        cbind(
                          chr = paste0(chr, "p"),
                          mut = mut_p,
                          ampl = ampl_p
                        ),
                        stringsAsFactors = FALSE
                      )
                    arm_ampl_mut <-
                      rbind(
                        arm_ampl_mut,
                        cbind(
                          chr = paste0(chr, "q"),
                          mut = mut_q,
                          ampl = ampl_q
                        ),
                        stringsAsFactors = FALSE
                      )
                  }
                }
                arm_ampl_mut[, 2:3] <-
                  apply(arm_ampl_mut[, 2:3], 2, as.numeric)
                (corP_arm <-
                    cor.test(
                      arm_ampl_mut$ampl,
                      log10(arm_ampl_mut$mut),
                      method = "pearson"
                    ))
                (corS_arm <-
                    cor.test(
                      arm_ampl_mut$ampl,
                      log10(arm_ampl_mut$mut),
                      method = "spearman"
                    ))
                
                parameters_arm <-
                  rbind(
                    parameters_arm,
                    c(
                      tumor_type = tumor_type,
                      type = "arm",
                      condition = condition,
                      corP = corP_arm$estimate,
                      p.valP = corP_arm$p.value,
                      corS = corS_arm$estimate,
                      p.valS = corS_arm$p.value
                    ),
                    stringsAsFactors = FALSE
                  )
                
                write_tsv(arm_ampl_mut,
                          file = paste0(results_table_path, tumor_type, "_chrArm.tsv"))
              }
              
              # produce plots
              if (condition != "deletions") {
                x <- log10(tier[tier$mutations_norm != 0, ]$mutations_norm)
                y <- tier[tier$mutations_norm != 0, ]$cna_freq_ampl
                
                corP <- try(cor.test(x, y, method = "pearson"))
                corS <- try(cor.test(x, y, method = "spearman"))
                
                parameters <- rbind(
                  parameters,
                  c(
                    tumor_type,
                    segment_length,
                    paste0(condition),
                    corP$estimate,
                    p.corP = corP$p.value,
                    corS$estimate,
                    p.corS = corS$p.value,
                    segment_cutoff_muts = segment_cutoff_muts,
                    segment_cutoff_CN = segment_cutoff
                  ),
                  stringsAsFactors = FALSE
                )
                (
                  p1 <- ggplot(tier[tier$mutations_norm != 0, ], aes(x = x, y = y)) +
                    geom_point() +
                    theme_classic() +
                    geom_smooth(method = "lm", alpha = 0.15) +
                    # geom_smooth(method = "gam", formula = y ~ s(x), alpha = 0.05, se = FALSE, size = 1, color = "black") +
                    ggtitle(
                      paste(tumor_type, "-", segment_length, "-", condition),
                      subtitle = paste(
                        "Pearson's R =",
                        signif(corP$estimate, digits = 3),
                        "\n",
                        "p-value =",
                        signif(corP$p.value, digits = 3),
                        "\nSpearman's rho =",
                        signif(corS$estimate, digits = 3),
                        "\n",
                        "p-value =",
                        signif(corS$p.value, digits = 3)
                      )
                    ) +
                    xlab("Log of diploid mutations") +
                    ylab("Amplification frequency")
                )
                
                system(
                  paste0(
                    "mkdir -p results/plots/02_tumorCorrelations/",
                    segment_length,
                    "Mbp_",
                    condition,
                    "/"
                  )
                )
                pdf(
                  file = paste0(
                    "results/plots/02_tumorCorrelations/",
                    segment_length,
                    "Mbp_",
                    condition,
                    "/",
                    tumor_type,
                    "_",
                    condition,
                    ".pdf"
                  )
                )
                print(p1)
                dev.off()
                
                if(condition == "non_expressed_GTEX" | 
                   condition == "expressed_GTEX" |
                   condition == "expressed_no0" |
                   condition == "non_expressed_no0" | 
                   condition == "haploinsufficient" |
                   condition == "non_haploinsufficient" |
                   condition == "removeWGD" |
                   condition == "all_mutations_mutsWithinCNA"){
                  x <- log10(tier[tier$mutations_norm != 0, ]$mutations_norm)
                  y <- tier[tier$mutations_norm != 0, ]$cna_freq_del
                  
                  corP <- try(cor.test(x, y, method = "pearson"))
                  corS <- try(cor.test(x, y, method = "spearman"))
                  
                  parameters <- rbind(
                    parameters,
                    c(
                      tumor_type,
                      segment_length,
                      paste0(condition, "_deletions"),
                      corP$estimate,
                      p.corP = corP$p.value,
                      corS$estimate,
                      p.corS = corS$p.value,
                      segment_cutoff_muts = segment_cutoff_muts,
                      segment_cutoff_CN = segment_cutoff
                    ),
                    stringsAsFactors = FALSE
                  )
                }
              }
              
              if (condition == "deletions") {
                x <- log10(tier[tier$mutations_norm != 0, ]$mutations_norm)
                y <- tier[tier$mutations_norm != 0, ]$cna_freq_del
                
                corP <- try(cor.test(x, y, method = "pearson"))
                corS <- try(cor.test(x, y, method = "spearman"))
                
                parameters <- rbind(
                  parameters,
                  c(
                    tumor_type,
                    segment_length,
                    paste0(condition),
                    corP$estimate,
                    p.corP = corP$p.value,
                    corS$estimate,
                    p.corS = corS$p.value,
                    segment_cutoff_muts = segment_cutoff_muts,
                    segment_cutoff_CN = segment_cutoff
                  ),
                  stringsAsFactors = FALSE
                )
                (
                  p1 <- ggplot(tier[tier$mutations_norm != 0, ], aes(x = x, y = y)) +
                    geom_point() +
                    theme_classic() +
                    geom_smooth(method = "lm", alpha = 0.15) +
                    # geom_smooth(method = "gam", formula = y ~ s(x), alpha = 0.05, se = FALSE, size = 1, color = "black") +
                    ggtitle(
                      paste(tumor_type, "-", segment_length, "-", condition),
                      subtitle = paste(
                        "Pearson's R =",
                        signif(corP$estimate, digits = 3),
                        "\n",
                        "p-value =",
                        signif(corP$p.value, digits = 3),
                        "\nSpearman's rho =",
                        signif(corS$estimate, digits = 3),
                        "\n",
                        "p-value =",
                        signif(corS$p.value, digits = 3)
                      )
                    ) +
                    xlab("Log of diploid mutations") +
                    ylab("Deletion frequency")
                )
                
                system(
                  paste0(
                    "mkdir -p results/plots/02_tumorCorrelations/",
                    segment_length,
                    "Mbp_",
                    condition,
                    "/"
                  )
                )
                pdf(
                  file = paste0(
                    "results/plots/02_tumorCorrelations/",
                    segment_length,
                    "Mbp_",
                    condition,
                    "/",
                    tumor_type,
                    "_",
                    condition,
                    ".pdf"
                  )
                )
                print(p1)
                dev.off()
              }
            }
            
            
          })
        }
        if (chr_arm == "yes") {
          colnames(parameters_chr) <- c("tumor_type",
                                        "type",
                                        "condition",
                                        "corP",
                                        "p.corP",
                                        "corS",
                                        "p.corS")
          colnames(parameters_arm) <- c("tumor_type",
                                        "type",
                                        "condition",
                                        "corP",
                                        "p.corP",
                                        "corS",
                                        "p.corS")
          parameters_chr[, 4:7] <-
            apply(parameters_chr[, 4:7], 2, as.numeric)
          parameters_arm[, 4:7] <-
            apply(parameters_arm[, 4:7], 2, as.numeric)
          
          write.table(
            parameters_chr,
            file = paste0(
              results_table_path,"02d_statistics_chromosome.txt"
            ),
            quote = F,
            row.names = F
          )
          write.table(
            parameters_arm,
            file = paste0(
              results_table_path,"02c_statistics_arm.txt"
            ),
            quote = F,
            row.names = F
          )
        }
      }
      
      colnames(parameters) <-
        c("tumorType",
          "segment_length",
          "condition",
          "corP",
          "p.corP",
          "corS",
          "p.corS",
          "cutoffMuts",
          "cutoffCN")
      parameters[, 4:7] <- apply(parameters[, 4:7], 2, as.numeric)
      
      if (length(segment_lengths) == 50) {
        write.table(
          parameters,
          file = paste0(
            results_table_path,"02b_statistics_from1to50Mbp_wConditions.txt"
          ),
          quote = F,
          row.names = F
        )
      }
    }
  }
}

# statistics ----
if (produce_statistics) {
  cat("\n >> (2) analyis: produce figures \n")
  
  
  fig1b_toPlot <- list()
  cat(" \n > Producing Fig. 1b \n\n")
  for(segment_cutoff in segment_cutoffs){
    print(segment_cutoff)
    
    toPlot <- data.frame()
    for (tumor_type in tumor_types) {
      print(tumor_type)
    scna <-
      read.delim(
        paste0(
          "data/FireBrowse_CNAs/",
          tumor_type,
          "_hg19_FireBrowse_totalCopyNumber_classifiedALL.tsv"
        ),
        header = T
      )
    scna_ampl <-
      scna[scna$Segment_Mean >= as.numeric(segment_cutoff)/100 & scna$length >= 500000, ]
    scna_del <-
      scna[scna$Segment_Mean <= -as.numeric(segment_cutoff)/100 & scna$length >= 500000, ]
    length_total_ampl <-
      sum(scna_ampl$length) / length(levels(factor(scna$patient_id)))
    length_total_del <-
      sum(scna_del$length) / length(levels(factor(scna$patient_id)))
    
    tier <-
      read.csv(
        paste0(
          "results/tables/02_produceStatistics/amplifications_",
          tumor_type,
          "_1Mbp_table",
          ifelse(segment_cutoff=="20","",paste0("_0",segment_cutoff)),".tsv"
        ),
        sep = "\t"
      )
    toPlot <- rbind(
      toPlot,
      cbind(
        tumor_type,
        cna_mean_freq_ampl = mean(tier$cna_freq_ampl),
        cna_mean_freq_del = mean(tier$cna_freq_del),
        cna_length_ampl = length_total_ampl,
        cna_length_del = length_total_del,
        mut_norm = log10(mean(tier$mutations_norm)),
        mut_mean_raw = mean(tier$mutations_raw),
        mut_tot_normPts = sum(tier$mutations_raw) /
          length(levels(factor(
            scna$patient_id
          ))),
        mut_tot = sum(tier$mutations_raw),
        patients = length(levels(factor(
          scna$patient_id
        )))
      )
    )
    }
    fig1b_toPlot[[segment_cutoff]] <- toPlot
  }
  
  for(segment_cutoff in segment_cutoffs){
    toPlot <- fig1b_toPlot[[segment_cutoff]]
    toPlot[, -1] <-  apply(toPlot[, -1], 2, as.numeric)
    
    # Figure 1b ----
    # correlation between mean mutation score and amplification frequency
    x <- as.numeric(toPlot$mut_norm)
    y <- as.numeric(toPlot$cna_mean_freq_ampl)
    
    corP <- cor.test(x, y, method = "pearson")
    corS <- cor.test(x, y, method = "spearman")
    
    (
      p1 <- ggplot(toPlot, aes(x = x, y = y)) +
        geom_point() +
        geom_text(aes(label = tumor_type), hjust = 0.5, vjust = 1.5) +
        theme_classic() +
        geom_smooth(method = "lm", alpha = 0.15) +
        ggtitle(
          paste("PANCANCER -", segment_cutoff),
          subtitle = paste(
            "Pearson's R =",
            signif(corP$estimate, digits = 3),
            "\n",
            "p-value =",
            signif(corP$p.value, digits = 3),
            "\nSpearman's rho =",
            signif(corS$estimate, digits = 3),
            "\n",
            "p-value =",
            signif(corS$p.value, digits = 3)
          )
        ) +
        xlab("Mean mutation score") +
        ylab("Mean amplification frequency")
    )
    
    pdf(
      paste0("results/plots/000_paper_plots/00_Fig1b_",segment_cutoff,".pdf"),
      width = 14,
      height = 8
    )
    print(p1)
    dev.off()
    
    png(paste0("results/plots/000_paper_plots/00_Fig1b_",segment_cutoff,".png"))
    print(p1)
    dev.off()
  }
  
  
  
  # correlating amplifications with deletions ----
  x <- as.numeric(toPlot$cna_mean_freq_ampl)
  y <- as.numeric(toPlot$cna_mean_freq_del)
  
  corP <- cor.test(x, y, method = "pearson")
  corS <- cor.test(x, y, method = "spearman")
  
  (
    p1 <- ggplot(toPlot, aes(x = x, y = y)) +
      geom_point() +
      geom_text(aes(label = tumor_type), hjust = 0.5, vjust = 1.5) +
      theme_classic() +
      geom_smooth(method = "lm", alpha = 0.15) +
      ggtitle(
        paste("PANCANCER"),
        subtitle = paste(
          "Pearson's R =",
          signif(corP$estimate, digits = 3),
          "\n",
          "p-value =",
          signif(corP$p.value, digits = 3),
          "\nSpearman's rho =",
          signif(corS$estimate, digits = 3),
          "\n",
          "p-value =",
          signif(corS$p.value, digits = 3)
        )
      ) +
      xlab("Mean amplification score") +
      ylab("Mean deletion frequency")
  )
  
  
  
  
  cat(" \n > Producing Fig. S2a \n\n")
  # Figure S1a ----
  # correlation between mean mutation score and deletion frequency
  x <- as.numeric(toPlot$mut_norm)
  y <- as.numeric(toPlot$cna_mean_freq_del)
  
  corP <- cor.test(x, y, method = "pearson")
  corS <- cor.test(x, y, method = "spearman")
  
  (
    p1 <- ggplot(toPlot, aes(x = x, y = y)) +
      geom_point() +
      geom_text(aes(label = tumor_type), hjust = 0.5, vjust = 1.5) +
      theme_classic() +
      geom_smooth(method = "lm", alpha = 0.15) +
      ggtitle(
        paste("PANCANCER"),
        subtitle = paste(
          "Pearson's R =",
          signif(corP$estimate, digits = 3),
          "\n",
          "p-value =",
          signif(corP$p.value, digits = 3),
          "\nSpearman's rho =",
          signif(corS$estimate, digits = 3),
          "\n",
          "p-value =",
          signif(corS$p.value, digits = 3)
        )
      ) +
      xlab("Mean mutation score") +
      ylab("Mean deletion frequency")
  )
  
  pdf(
    "results/plots/000_paper_plots/00_SupplementaryFig2a.pdf",
    width = 14,
    height = 8
  )
  print(p1)
  dev.off()
  
  cat(" \n > Producing Fig. 1f \n\n")
  # Figure 1f ----
  # correlation between mutation burden and correlation estimate
  corParameters <-
    read.table(
      "results/tables/02_produceStatistics/02b_statistics_from1to50Mbp_wConditions.txt",
      header = T
    )
  corParameters <-
    corParameters[corParameters$segment_length == 36 &
                    corParameters$condition == "amplifications", ]
  colnames(corParameters)[1] <- "tumor_type"
  
  join <- left_join(corParameters, toPlot, by = "tumor_type")
  
  x <- log10(as.numeric(join$mut_mean_raw))
  y <- as.numeric(join$corS)
  
  corP <- cor.test(x, y, method = "pearson")
  corS <- cor.test(x, y, method = "spearman")
  
  (
    p1 <- ggplot(join, aes(x = x, y = y)) +
      geom_point() +
      geom_text(aes(label = tumor_type), hjust = 0.5, vjust = 1.5) +
      theme_classic() +
      geom_smooth(method = "lm", alpha = 0.15) +
      ggtitle(
        paste("PANCANCER"),
        subtitle = paste(
          "Pearson's R =",
          signif(corP$estimate, digits = 3),
          "\n",
          "p-value =",
          signif(corP$p.value, digits = 3),
          "\nSpearman's rho =",
          signif(corS$estimate, digits = 3),
          "\n",
          "p-value =",
          signif(corS$p.value, digits = 3)
        )
      ) +
      xlab("Mutation burden (log10)") +
      ylab("Spearman Correlation")
  )
  
  pdf(
    "results/plots/000_paper_plots/00_Fig1f.pdf",
    width = 14,
    height = 8
  )
  print(p1)
  dev.off()
  
  cat(" \n > Producing Fig. S1b \n\n")
  # Figure S1b ----
  # correlation between amplification burden and correlation estimate
  x <- as.numeric(join$cna_mean_freq_ampl)
  y <- as.numeric(join$corS)
  
  corP <- cor.test(x, y, method = "pearson")
  corS <- cor.test(x, y, method = "spearman")
  
  (
    p1 <- ggplot(join, aes(x = x, y = y)) +
      geom_point() +
      geom_text(aes(label = tumor_type), hjust = 0.5, vjust = 1.5) +
      theme_classic() +
      geom_smooth(method = "lm", alpha = 0.15) +
      ggtitle(
        paste("PANCANCER"),
        subtitle = paste(
          "Pearson's R =",
          signif(corP$estimate, digits = 3),
          "\n",
          "p-value =",
          signif(corP$p.value, digits = 3),
          "\nSpearman's rho =",
          signif(corS$estimate, digits = 3),
          "\n",
          "p-value =",
          signif(corS$p.value, digits = 3)
        )
      ) +
      xlab("Amplification Burden (log10)") +
      ylab("Spearman Correlation")
  )
  
  pdf(
    "results/plots/000_paper_plots/00_SupplementaryFig2b.pdf",
    width = 14,
    height = 8
  )
  print(p1)
  dev.off()
  
  cat(" \n > Producing Fig. 1e \n\n")
  # Figure 1e ----
  # curve of correlation estimates varying correlation length
  levelChr <-
    read.delim("results/tables/02_produceStatistics/02d_statistics_chromosome.txt",
               sep = " ")
  levelChr <-
    cbind(levelChr, segment_length = rep("70", dim(levelChr)[1]))
  
  levelArm <-
    read.delim("results/tables/02_produceStatistics/02c_statistics_arm.txt",
               sep = " ")
  levelArm <-
    cbind(levelArm, segment_length = rep("60", dim(levelArm)[1]))
  
  levelGene  <-
    read.delim("results/tables/02_produceStatistics/02a_statistics_gene.txt",
               sep = " ")
  levelGene <-
    cbind(levelGene, segment_length = rep("0", dim(levelGene)[1]))
  
  final <- rbind(levelArm, levelChr,
                 levelGene)
  
  statistics_acd <- data.frame()
  for (i in levels(factor(final$segment_length))) {
    statistics_acd <- rbind(statistics_acd,
                            cbind(
                              i,
                              corP = mean(as.numeric(final[final$segment_length == i, ]$corP)),
                              p.corP = mean(as.numeric(final[final$segment_length == i, ]$p.corP)),
                              corS = mean(as.numeric(final[final$segment_length == i, ]$corS)),
                              p.corS = mean(as.numeric(final[final$segment_length == i, ]$p.corS))
                            ))
  }
  
  corParameters <-
    read.table(
      "results/tables/02_produceStatistics/02b_statistics_from1to50Mbp_wConditions.txt",
      header = T
    )
  corParameters <-
    corParameters[corParameters$condition == "amplifications", ]
  
  table_selected_tumor_types <-
    table(corParameters$tumorType, corParameters$p.corS <= 0.05)
  
  # selected_tumor_types <- names(table_selected_tumor_types[,2])[table_selected_tumor_types[,2] >= 0] # uncomment this for an unbiased analysis (all tumor types)
  selected_tumor_types <-
    names(table_selected_tumor_types[, 2])[table_selected_tumor_types[, 2] >= 30] # include only highly significant cancer types
  
  
  corParameters <-
    corParameters[corParameters$tumorType %in% selected_tumor_types, ]
  final <- final[final$tumor_type %in% selected_tumor_types, ]
  
  statistics_b <- data.frame()
  for (i in levels(factor(corParameters$segment_length))) {
    statistics_b <- rbind(statistics_b,
                          cbind(
                            i,
                            corP = mean(corParameters[corParameters$segment_length == i, ]$corP),
                            p.corP = mean(corParameters[corParameters$segment_length == i, ]$p.corP),
                            corS = mean(corParameters[corParameters$segment_length == i, ]$corS),
                            p.corS = mean(corParameters[corParameters$segment_length == i, ]$p.corS)
                          ))
  }
  
  statistics <- rbind(statistics_b, statistics_acd)
  
  # statistics <- left_join(data.frame(i = as.factor(c(0,1,2,4,6,8,10,12,14,16,18,20,22,24,28,30,32,34,36,40,44,50,60,70))), statistics) # optional
  statistics <- as.data.frame(apply(statistics, 2, as.numeric))
  
  p1 <- ggplot(statistics, aes(y = corS, x = -i)) +
    geom_point() +
    stat_smooth(method = "lm",
                formula = y ~ x + I(x ^ 2),
                size = 1) +
    geom_vline(
      xintercept = -36,
      linetype = "dashed",
      color = "red",
      size = 0.5
    ) +
    ggtitle("PANCANCER - Spearman's estimate across different chromosomes lengths") +
    theme_classic()
  
  pdf(
    "results/plots/000_paper_plots/00_Fig1e.pdf",
    width = 14,
    height = 8
  )
  print(p1)
  dev.off()
  
  cat(" \n > Producing Fig. 1c \n\n")
  # Figure 1c ----
  corParameters <-
    read.table(
      "results/tables/02_produceStatistics/02b_statistics_from1to50Mbp_wConditions.txt",
      header = T
    )
  corParameters <-
    corParameters[corParameters$segment_length == 36 &
                    corParameters$condition == "amplifications", ]
  corParameters$tumorType <- factor(corParameters$tumorType,
                                    levels = c(corParameters[order(corParameters[corParameters$condition == "amplifications", ]$corS, decreasing = T), ]$tumorType))
  
  amplifications <-
    corParameters[corParameters$condition == "amplifications", ]
  p1 <-
    ggplot(amplifications, aes(
      x = as.factor(tumorType),
      y = corS,
      fill = as.factor(condition)
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Blues") +
    theme(legend.position = "none")
  p2 <-
    ggplot(amplifications, aes(
      x = as.factor(tumorType),
      y = log10(p.corS),
      fill = condition
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Reds") +
    geom_hline(yintercept = log10(0.05)) +
    theme(legend.position = "none")
  
  pdf(
    "results/plots/000_paper_plots/00_Fig1c.pdf",
    width = 14,
    height = 8
  )
  print(ggarrange(p1, p2, nrow = 2, ncol = 1))
  dev.off()
  
  # filter only positive and significant cancer type (Spearman's p-value <= 0.05)
  corParameters <-
    read.table(
      "results/tables/02_produceStatistics/02b_statistics_from1to50Mbp_wConditions.txt",
      header = T
    )
  corParameters <-
    corParameters[corParameters$tumorType %in% na.omit(corParameters[corParameters$condition == "amplifications" &
                                                                       corParameters$p.corS <= 0.05 &
                                                                       corParameters$segment_length == 36, ]$tumorType), ]
  corParameters <-
    corParameters[corParameters$segment_length == 36, ]
  
  # estimate the differences between groups
  wilcox.test(
    corParameters[corParameters$condition == "amplifications", ]$corS,
    corParameters[corParameters$condition == "deletions", ]$corS,
    paired = T,
    alternative = "greater"
  )
  wilcox.test(
    corParameters[corParameters$condition == "expressed_no0", ]$corS,
    corParameters[corParameters$condition == "non_expressed_no0", ]$corS,
    paired = T,
    alternative = "greater"
  )
  wilcox.test(
    corParameters[corParameters$condition == "CADD_moderatelyDamaging_phred", ]$corS,
    corParameters[corParameters$condition == "CADD_highlyDamaging_phred", ]$corS,
    paired = T,
    alternative = "less"
  )
  wilcox.test(
    corParameters[corParameters$condition == "polyphen_moderatelyDamaging", ]$corS,
    corParameters[corParameters$condition == "polyphen_highlyDamaging", ]$corS,
    paired = T,
    alternative = "less"
  )
  wilcox.test(
    corParameters[corParameters$condition == "aggregation_causing", ]$corS,
    corParameters[corParameters$condition == "non_aggregation_causing", ]$corS,
    paired = T,
    alternative = "less"
  )
  wilcox.test(
    corParameters[corParameters$condition == "haploinsufficient", ]$corS,
    corParameters[corParameters$condition == "non_haploinsufficient", ]$corS,
    paired = T,
    alternative = "greater"
  )
  wilcox.test(
    corParameters[corParameters$condition == "haploinsufficient_GHIS", ]$corS,
    corParameters[corParameters$condition == "non_haploinsufficient_GHIS", ]$corS,
    paired = T,
    alternative = "greater"
  )
  wilcox.test(corParameters[corParameters$condition == "amplifications", ]$corS, corParameters[corParameters$condition == "remove_OG", ]$corS)
  wilcox.test(corParameters[corParameters$condition == "amplifications", ]$corS, corParameters[corParameters$condition == "remove_TSG", ]$corS)
  wilcox.test(corParameters[corParameters$condition == "amplifications", ]$corS, corParameters[corParameters$condition == "remove_BOTH", ]$corS)
  
  
  amplifications <-
    corParameters[corParameters$condition == "amplifications", ]
  p1 <-
    ggplot(amplifications, aes(
      x = as.factor(tumorType),
      y = corS,
      fill = as.factor(condition)
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Blues")
  p2 <-
    ggplot(amplifications, aes(
      x = as.factor(tumorType),
      y = log10(p.corS),
      fill = condition
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Reds") +
    geom_hline(yintercept = log10(0.05))
  
  print(ggarrange(p1, p2, nrow = 2, ncol = 1))
  
  cat(" \n > Producing Fig. 1d \n\n")
  # Figure 1d ----
  deletions <-
    corParameters[corParameters$condition == "amplifications" |
                    corParameters$condition == "deletions", ]
  p1 <-
    ggplot(deletions, aes(
      x = as.factor(tumorType),
      y = corS,
      fill = as.factor(condition)
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Blues")
  p2 <-
    ggplot(deletions, aes(
      x = as.factor(tumorType),
      y = log10(p.corS),
      fill = condition
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Reds")
  
  pdf(
    "results/plots/000_paper_plots/00_Fig1d.pdf",
    width = 14,
    height = 8
  )
  print(
    ggplot(deletions, aes(
      x = corS,
      y = condition,
      fill = as.factor(condition)
    )) +
      geom_boxplot() +
      theme_classic() +
      scale_fill_brewer(palette = "Reds") +
      geom_jitter()
  )
  print(ggarrange(p1, p2, nrow = 2, ncol = 1))
  dev.off()
  
  
  
  cat(" \n > Producing Fig. 2a-b \n\n")
  # figure 2a-b ----
  haploinsuf <-
    corParameters[corParameters$condition == "amplifications" |
                    corParameters$condition == "haploinsufficient" |
                    corParameters$condition == "non_haploinsufficient" |
                    corParameters$condition == "haploinsufficient_GHIS" |
                    corParameters$condition == "non_haploinsufficient_GHIS", ]
  p1 <-
    ggplot(haploinsuf, aes(
      x = as.factor(tumorType),
      y = corS,
      fill = as.factor(condition)
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Blues")
  p2 <-
    ggplot(haploinsuf, aes(
      x = as.factor(tumorType),
      y = log10(p.corS),
      fill = condition
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Reds")
  
  pdf(
    "results/plots/000_paper_plots/00_Fig2a.pdf",
    width = 10,
    height = 8
  )
  print(ggarrange(p1, p2, nrow = 2, ncol = 1))
  dev.off()
  
  cat(" \n > Producing Fig. 2c-d \n\n")
  # figure 2c-d----
  haploinsuf <-
    corParameters[corParameters$condition == "amplifications" |
                    corParameters$condition == "expressed_no0" |
                    corParameters$condition == "non_expressed_no0", ]
  p1 <-
    ggplot(haploinsuf, aes(
      x = as.factor(tumorType),
      y = corS,
      fill = as.factor(condition)
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Blues")
  p2 <-
    ggplot(haploinsuf, aes(
      x = as.factor(tumorType),
      y = log10(p.corS),
      fill = condition
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Reds")
  
  pdf(
    "results/plots/000_paper_plots/00_Fig2c.pdf",
    width = 10,
    height = 8
  )
  print(ggarrange(p1, p2, nrow = 2, ncol = 1))
  dev.off()
  
  
  
  cat(" \n > Producing Fig. 2e-f \n\n")
  # figure 2e-f ----
  pdf(
    "results/plots/000_paper_plots/00_Fig2e.pdf",
    width = 10,
    height = 8
  )
  
  polycadd <-
    corParameters[corParameters$condition == "amplifications" |
                    corParameters$condition == "polyphen_highlyDamaging" |
                    corParameters$condition == "polyphen_moderatelyDamaging"  |
                    corParameters$condition == "CADD_moderatelyDamaging" |
                    corParameters$condition == "CADD_highlyDamaging", ]
  p1 <-
    ggplot(polycadd, aes(
      x = as.factor(tumorType),
      y = corS,
      fill = condition
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Blues")
  p2 <-
    ggplot(polycadd, aes(
      x = as.factor(tumorType),
      y = log10(p.corS),
      fill = condition
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Reds")
  
  print(ggarrange(p1, p2, nrow = 2, ncol = 1))
  
  cadd <-
    corParameters[corParameters$condition == "amplifications" |
                    corParameters$condition == "CADD_moderatelyDamaging" |
                    corParameters$condition == "CADD_highlyDamaging", ]
  p1 <-
    ggplot(cadd, aes(
      x = as.factor(tumorType),
      y = corS,
      fill = condition
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Blues")
  p2 <-
    ggplot(cadd, aes(
      x = as.factor(tumorType),
      y = log10(p.corS),
      fill = condition
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Reds")
  
  print(ggarrange(p1, p2, nrow = 2, ncol = 1))
  
  polyphen <-
    corParameters[corParameters$condition == "amplifications" |
                    corParameters$condition == "polyphen_highlyDamaging" |
                    corParameters$condition == "polyphen_moderatelyDamaging", ]
  p1 <-
    ggplot(polyphen, aes(
      x = as.factor(tumorType),
      y = corS,
      fill = condition
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Blues")
  p2 <-
    ggplot(polyphen, aes(
      x = as.factor(tumorType),
      y = log10(p.corS),
      fill = condition
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Reds")
  
  print(ggarrange(p1, p2, nrow = 2, ncol = 1))
  dev.off()
  
  
  
  cat(" \n > Producing Fig. 2g-h \n\n")
  # figure 2g-h----
  haploinsuf <-
    corParameters[corParameters$condition == "amplifications" |
                    corParameters$condition == "synonymous" |
                    corParameters$condition == "non_synonymous", ]
  p1 <-
    ggplot(haploinsuf, aes(
      x = as.factor(tumorType),
      y = corS,
      fill = as.factor(condition)
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Blues")
  p2 <-
    ggplot(haploinsuf, aes(
      x = as.factor(tumorType),
      y = log10(p.corS),
      fill = condition
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Reds")
  
  pdf(
    "results/plots/000_paper_plots/00_Fig2g.pdf",
    width = 10,
    height = 8
  )
  print(ggarrange(p1, p2, nrow = 2, ncol = 1))
  dev.off()
  
  
  
  cat(" \n > Producing Fig. 2i-j \n\n")
  # figure 2i-j ----
  aggregat <-
    corParameters[corParameters$condition == "amplifications" |
                    corParameters$condition == "aggregation_causing" |
                    corParameters$condition == "non_aggregation_causing", ]
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
  
  pdf(
    "results/plots/000_paper_plots/00_Fig2i.pdf",
    width = 10,
    height = 8
  )
  print(ggarrange(p1, p2, nrow = 2, ncol = 1))
  dev.off()
  
  
  
  cat(" \n > Producing Fig. 2k-l \n\n")
  # figure 2k-l ----
  cancergenes <-
    corParameters[corParameters$condition == "amplifications" |
                    corParameters$condition == "remove_OG" |
                    corParameters$condition == "remove_TSG" |
                    corParameters$condition == "remove_BOTH", ]
  p1 <-
    ggplot(cancergenes, aes(
      x = as.factor(tumorType),
      y = corS,
      fill = condition
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Blues")
  p2 <-
    ggplot(cancergenes, aes(
      x = as.factor(tumorType),
      y = log10(p.corS),
      fill = condition
    )) +
    geom_bar(color = "black",
             stat = "identity",
             position = position_dodge()) +
    theme_classic() + scale_fill_brewer(palette = "Reds")
  
  pdf(
    "results/plots/000_paper_plots/00_Fig2k.pdf",
    width = 10,
    height = 8
  )
  print(ggarrange(p1, p2, nrow = 2, ncol = 1))
  dev.off()
}

cat(
  "\n\n OUTPUT of the script: \n \t (1) raw tables path: results/plots/02_tumorCorrelations/\n"
)
cat("\t (2) correlation tables path: results/tables/02_produceStatistics/\n")
cat("\t (3) paper-like plots path: results/plots/000_paper_plots/ \n\n")

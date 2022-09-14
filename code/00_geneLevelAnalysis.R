# GENE-LEVEL ANALYSIS
# This script integrates data from TCGA to compute gene specific mutation score
# and amplification frequency

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
  library(foreach)
  library(doSNOW)
  library(snow)
  library(optparse)
})

option_list = list(
  make_option(
    c("-t", "--tables"),
    type = "character",
    default = "n",
    help = "Options are: (y/[n])
              ATTENTION:  If 'y', it may take several hours and requires parallelization 
              (see cores parameter within the Rscript and set according to your machine). 
              Set 'n' to skip this step, but only if you already downloaded the preprocessed 
              data contained results.zip folder.",
    metavar = ""
  ),
  make_option(
    c("-s", "--statistics"),
    type = "character",
    default = "y",
    help = "Options are: ([y]/n)
              If 'y' calculates gene level correlations; otherwise it skips this part.",
    metavar = ""
  )
)


opt_parser = OptionParser(option_list = option_list)

opt = parse_args(opt_parser)


if (opt$tables == "n" & opt$statistics == "y") {
  cat("\n\n >> You chose default options: \n\t(1) --tables 'n';\n\t(2) --statistics 'y' \n\n")
}

if (is.null(opt$tables) | is.null(opt$statistics)) {
  print_help(opt_parser)
  stop("please specify the analysis you want to perform!", call. = FALSE)
}

if (!any(opt$tables %in% c("y", "n"))) {
  print_help(opt_parser)
  stop("typo in the analysis flag, plase see above for the available options!",
       call. = FALSE)
}

cat(
  "\n\n > This script \n\t(1) estimates gene amplification frequency and mu score (takes several hours and cores);\n\t(2) produces gene-level correlations \n\n\n"
)

produce_tables <- F
if (opt$tables == "y") {
  produce_tables <- T
}
produce_statistics <- T
if (opt$statistics == "n") {
  produce_statistics <- F
}

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
  "UCS",
  "SKCM"
)

setwd("../")

# path for gene-level analysis output
results_table_path <- "results/tables/00_geneLevelAnalysis/"
system(paste0("mkdir -p ", results_table_path))

results_plots_path <- "results/plots/00_geneLevelAnalysis/"
system(paste0("mkdir -p ", results_plots_path))

# WARNING (only for (1) analysis)
cores <- 22 # set cores based on your core availability
# only used when produce_table == T

if (produce_tables) {
  cat("\n >> (1) analyis: estimates gene amplification frequency and mu score \n")
  # it generates a table for each cancer type
  for (tumor_type in tumor_types) {
    cat("\n\n >>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<< \n")
    cat(">      ----------", paste0(tumor_type, " ---------- \n"))
    
    # load SNVs file
    cat("> Load Mutation File \n")
    snv <-
      read.csv(file = paste0("data/FireBrowse_SNVs/", tumor_type, "_mutations.csv"))
    snv <-
      snv[!duplicated(snv[, c(6, 7, 8, 9, 17)]),] # remove dupicated mutations in the same patient
    cat("  done!\n\n")
    
    # load CNAs file
    cat("> Load SCNA File \n")
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
    scna <- scna[scna$patient_id %in% snv$patient_id,]
    common_patients <-
      levels(as.factor(snv[snv$patient_id %in% scna$patient_id,]$patient_id))
    scna <-
      scna[scna$Segment_Mean >= 0.2 | scna$Segment_Mean <= -0.2,]
    cat("  done!\n\n")
    cat("> common patients:", length(common_patients), "\n")
    
    proteins <-
      read.table("data/misc/protein-coding_gene.txt", header = FALSE)
    # downloaded from: ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/
    coding_regions <-
      read.table("data/misc/CCDS.current_hg19.txt",
                 header = TRUE,
                 sep = "\t")
    # filter just protein coding genes
    coding_regions <-
      coding_regions[coding_regions$gene %in% proteins$V1, ]
    coding_regions <-
      coding_regions[coding_regions$ccds_status == "Public", ]
    coding_regions$length <-
      as.numeric(coding_regions$cds_to) - as.numeric(coding_regions$cds_from)
    coding_regions <-
      coding_regions[order(coding_regions$length, decreasing = T), ]
    coding_regions <-
      coding_regions[!duplicated(coding_regions$gene), ]
    total_pts <- length(levels(factor(scna$patient_id)))
    
    cat("> Load TPM File \n")
    tpm <-
      readRDS(file = paste0("data/TCGA_tpm/", tumor_type, "_tpm.rds.gz"))
    cat("  done!\n\n")
    
    # amplifications analysis
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    
    cat("> start parallalization with", cores, "cores \n")
    # for each gene assess the number of
    results <- foreach(chr = 1:22, .combine = "rbind") %dopar% {
      temp_coding <- coding_regions[coding_regions$chromosome == chr, ]
      temp_mut <- snv_filt[snv_filt$Chromosome == chr, ]
      temp_cna <- scna[scna$Chromosome == chr, ]
      n_genes <- dim(temp_coding)[1]
      result <- data.frame()
      for (i in 1:nrow(temp_coding)) {
        gene_name <- temp_coding$gene[i]
        # extract CNA in this gene
        cna_gene <-
          temp_cna[as.numeric(temp_cna$Start) <= as.numeric(temp_coding[i, ]$cds_from) &
                     as.numeric(temp_cna$End) >= as.numeric(temp_coding[i, ]$cds_to), ]
        # assess the frequency of the amplification
        amplified_pts <- levels(factor(cna_gene$patient_id))
        CN_pts <-
          common_patients[!common_patients %in% amplified_pts]
        
        gene_tpm <-
          tpm[tpm$patient_id %in% CN_pts &
                tpm$`Approved symbol` == gene_name, ]
        gene_tpm <- mean(gene_tpm$tpm_counts, na.rm = T)
        
        freq_ampl <-
          length(levels(factor(cna_gene$patient_id))) / as.numeric(total_pts)
        
        n_mut_all <-
          dim(temp_mut[as.numeric(temp_mut$Start_Position) >= as.numeric(temp_coding[i, ]$cds_from) &
                         as.numeric(temp_mut$Start_Position) <= as.numeric(temp_coding[i, ]$cds_to), ])[1] /
          total_pts
        gene_length <-
          as.numeric(temp_coding$cds_to[i]) - as.numeric(temp_coding$cds_from[i])
        n_mut_norm <- n_mut_all / gene_length
        
        n_pt_CN <- total_pts - length(amplified_pts)
        temp_mut_CN <-
          temp_mut[!temp_mut$patient_id %in% amplified_pts, ]
        mutations_CN <-
          temp_mut_CN[as.numeric(temp_mut_CN$Start_Position) >= as.numeric(temp_coding[i, ]$cds_from) &
                        as.numeric(temp_mut_CN$Start_Position) <= as.numeric(temp_coding[i, ]$cds_to), ]
        n_mut_all_CN <- dim(mutations_CN)[1] / n_pt_CN
        n_mut_CN_norm <- n_mut_all_CN / temp_coding[i, "length"]
        
        result <- rbind(
          result,
          cbind(
            Chromosome = chr,
            Gene = gene_name,
            Amplification = freq_ampl,
            CN_Pts = n_pt_CN,
            GeneTPM = gene_tpm,
            GeneLength = gene_length,
            MutationsCN = n_mut_all_CN,
            MutationsCN_norm = n_mut_CN_norm
          )
        )
      }
      return(result)
    }
    
    stopCluster(cl)
    
    # remove duplicate genes or NA (if present)
    results <-
      results[!duplicated(paste0(results$Gene, results$GeneLegth)),]
    results <-
      results[!(
        is.na(results$Amplification) |
          is.na(results$MutationsCN_norm) |
          is.na(results$GeneTPM)
      ),]
    
    # produce the tables
    write.csv(results,
              file = paste0(results_table_path, tumor_type, "_geneLevel_TRY.csv"))
  }
}

if (produce_statistics) {
  cat("\n >> (2) analyis: produces gene-level correlations \n")
  
  # produce correlation estimates at gene level ----
  statistics <- data.frame()
  
  for (tumor_type in tumor_types) {
    results <-
      read.csv(paste0(results_table_path, tumor_type, "_geneLevel.csv"))
    
    results <- results[results$MutationsDiploid_norm != 0, ]
    results <- results[results$Amplification != 0, ]
    results <-
      results[!duplicated(paste0(results$Gene, results$GeneLegth)), ]
    
    y <- as.numeric(results$Amplification)
    x <- log10(as.numeric(results$MutationsDiploid_norm))
    
    suppressWarnings({
      corP <- cor.test(x, y, method = "pearson")
      corS <- cor.test(x, y, method = "spearman")
    })
    
    statistics <- rbind(
      statistics,
      cbind(
        tumor_type,
        type = "gene",
        condition = "amplifications",
        corP = corP$estimate,
        corS = corS$estimate,
        p.corP = corP$p.value,
        p.corS = corS$p.value
      )
    )
    
    p1 <- ggplot(results, aes(y = y, x = x)) +
      geom_point() +
      theme_classic() +
      geom_smooth(method = "lm", alpha = 0.15) +
      xlab("Mutation score") +
      ylab("Amplification frequency") +
      ggtitle(
        paste0(tumor_type, " - All genes (", dim(results)[1], ")"),
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
      )
    pdf(file = paste0(results_plots_path, tumor_type, "_geneAnalysis.pdf"))
    print(p1)
    dev.off()
  }
  
  cat("\n\n > Spearman's correlations for GENE-LEVEL analysis: \n")
  print(statistics)[-c(3, 4, 6)]
  
  results_table_statistics_path <-
    "results/tables/02_produceStatistics/"
  system(paste0("mkdir -p ", results_table_statistics_path))
  write.table(
    statistics,
    file = paste0(results_table_statistics_path, "02a_statistics_gene.txt"),
    quote = F,
    row.names = F
  )
}

cat("\n\n > OUTPUT of the script: \n \t (1) raw tables path:",
    results_table_path,
    "\n")
cat("\t (2) correlation plots path:", results_plots_path, "\n")
cat("\t (3) gene-level correlations path:",
    results_table_statistics_path,
    "\n\n")
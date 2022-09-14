# TIMING ANALYSIS using CNAqc

suppressMessages({
  library(ggridges)
  library(tidyr)
  library(readr)
  library(dbplyr)
  library(devtools)
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(readr)
})


option_list = list(
  make_option(
    c("-r", "--runCNAqc"),
    type = "character",
    default = "n",
    
    help = "Options are: [y/(n)]
              set 'y' if you want to produce PCAWG
              ATTENTION: --runCNAqc set to 'y' requires lots of computational power and lot of time",
    metavar = ""
  ),
  make_option(
    c("-m", "--mock"),
    type = "character",
    default = "y",
    help = "Options are: [(y)/n]
              Performs the exact same analysis with a restricted pool of mutations (n = 100,000)",
    metavar = ""
  ),
  make_option(
    c("-c", "--clonal"),
    type = "character",
    default = "n",
    help = "Options are: [(y)/n]
              Performs the analysis only on clonal mutations",
    metavar = ""
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



if (opt$runCNAqc == "n" & opt$mock == "y") {
  cat("\n\n >> You chose default options: \n\t(1) --runCNAqc 'n';\n\t(2) --mock 'y' \n\n")
}

if (is.null(opt$runCNAqc) | is.null(opt$mock)) {
  print_help(opt_parser)
  stop("please specify the analysis you want to perform!", call. = FALSE)
}

cat(
  "\n\n > This script \n\t(1) estimates gene amplification frequency and mu score (takes several hours and cores);\n\t(2) produces gene-level correlations \n\n\n"
)

setwd("../")
system(paste0("mkdir -p ", getwd(), "/results/plots/04_timingAnalysis/"))

# set to TRUE if you want to re-run the CNAqc phasing analysis on the whole PCAWG dataset
# otherwise FALSE to use pre-computed PCAWG dataset and only run the timing analysis
# (it may take several hours)
run_CNAqc_analysis <- F
if (opt$runCNAqc == "y") {
  run_CNAqc_analysis <- T
}

# set to TRUE to perform the timing analysis on a mock dataset of 100,000 mutations
# (the actual dataset is about 5GB - may take several minutes to load into R)
test <- F
if (opt$mock == "y") {
  test <- T
}

# set to TRUE if you want to perform the analysis only on clonal mutations, otherwise
# all mutations will be included (clonal and subclonal)
clonal <- F
if (opt$clonal == "y") {
  clonal <- T
}


if (run_CNAqc_analysis & !test) {
  # install.packages("devtools")
  if (!require("CNAqc")) {
    devtools::install_github("caravagnalab/CNAqc")
  }
  library(CNAqc)
  
  # load the list of patients within PCAWG database
  patientsPCAWG <-
    read.table(file = "data/misc/PCAWG_n_snv_cna.txt", header = T)
  load("data/misc/CNAqc_functions.RData")
  
  # run the advanced_phasing() analysis (CNAqc) for each PCAWG patient
  PCAWG_PATH <-
    paste0("data/clonal_analysis_PCAWG/")
  system("mkdir -p results/tables/04_timingAnalysis/tmp")
  
  for (pcawg_id in patientsPCAWG$pt) {
    print(pcawg_id)
    pt <- readRDS(paste0(PCAWG_PATH, pcawg_id, "/fit.rds"))
    
    try({
      library(CNAqc)
      library(dplyr)
      library(ggplot2)
      library(readr)
      
      # run phasing analysis using CNAqc https://caravagnalab.github.io/CNAqc/
      pt <- advanced_phasing(pt)
      CCF_table <- pt$phasing
      
      CCF_table$purity <-
        patientsPCAWG[patientsPCAWG$pt == pcawg_id, ]$purity
      CCF_table$ploidy <-
        patientsPCAWG[patientsPCAWG$pt == pcawg_id, ]$ploidy
      CCF_table$tumor_type <-
        patientsPCAWG[patientsPCAWG$pt == pcawg_id, ]$tumor_type
      
      write_rds(
        CCF_table,
        paste0(
          "results/tables/04_timingAnalysis/tmp/",
          pcawg_id,
          "_CCF.rds"
        )
      )
    })
  }
  
  # merge single patient phased mutation dataset in a single dataset
  mutations_wCCF <- data.frame(stringsAsFactors = F)
  n <- 0
  for (file in paste0(patientsPCAWG$pt, "_CCF.rds")) {
    try({
      pt_file <-
        readRDS(paste0("results/tables/04_timingAnalysis/tmp/", file))
      # print(colnames(pt_file))
      mutations_wCCF <- rbind(mutations_wCCF, pt_file)
      n <- n + 1
    })
  }
  write.table(
    mutations_wCCF,
    file = paste0("data/misc/PCAWG_mutations_wCCF_allTumorTypes.tsv"),
    sep = "\t"
  )
}

if (test) {
  mutations_wCCF <-
    read_tsv(paste0("data/misc/PCAWG_mutations_wCCF_mockFile.tsv"))
} else{
  # ATTENTION: big dataset (5 GB) !
  mutations_wCCF <-
    read_tsv(paste0("data/misc/PCAWG_mutations_wCCF_allTumorTypes.tsv"))
}

mutations_wCCF$multiplicity <-
  as.factor(mutations_wCCF$multiplicity)

if (clonal) {
  # filter only clonal mutations (CCF >= 1)
  mutations_wCCF <- mutations_wCCF[mutations_wCCF$CCF >= 1, ]
}

# define coding and non-coding mutations ----
mutations_wCCF$multiplicity <-
  as.numeric(as.character(mutations_wCCF$multiplicity))

intergenic_region <-
  mutations_wCCF %>% filter(consequence_type == "intergenic_region")
intergenic_region_mixed <-
  mutations_wCCF %>% filter(
    consequence_type == "intergenic_region:upstream_gene_variant" |
      consequence_type == "intergenic_region:downstream_gene_variant" |
      consequence_type == "downstream_gene_variant:intergenic_region" |
      consequence_type == "upstream_gene_variant:intergenic_region" |
      consequence_type == "intergenic_region:upstream_gene_variant:downstream_gene_variant" |
      consequence_type == "intergenic_region:downstream_gene_variant:upstream_gene_variant" |
      consequence_type == "downstream_gene_variant:intergenic_region:upstream_gene_variant" |
      consequence_type == "upstream_gene_variant:intergenic_region:downstream_gene_variant" |
      consequence_type == "downstream_gene_variant:upstream_gene_variant:intergenic_region" |
      consequence_type == "upstream_gene_variant:downstream_gene_variant:intergenic_region"
  )

exon_variant  <-
  mutations_wCCF %>% filter(consequence_type == "exon_variant")
exon_variant_mixed  <-
  mutations_wCCF %>% filter(
    consequence_type == "intron_variant:exon_variant" |
      consequence_type == "exon_variant:intron_variant" |
      consequence_type == "downstream_gene_variant:exon_variant" |
      consequence_type == "upstream_gene_variant:exon_variant" |
      consequence_type == "intron_variant:downstream_gene_variant:exon_variant" |
      consequence_type == "exon_variant:upstream_gene_variant" |
      consequence_type == "exon_variant:downstream_gene_variant" |
      consequence_type == "intron_variant:upstream_gene_variant:exon_variant" |
      consequence_type == "downstream_gene_variant:intron_variant:exon_variant" |
      consequence_type == "intron_variant:exon_variant:downstream_gene_variant" |
      consequence_type == "intron_variant:exon_variant:upstream_gene_variant"
  )

missense_variant <-
  mutations_wCCF %>% filter(consequence_type == "missense_variant")
missense_variant_mixed <-
  mutations_wCCF %>% filter(
    consequence_type == "missense_variant:exon_variant" |
      consequence_type == "missense_variant:downstream_gene_variant" |
      consequence_type == "intron_variant:missense_variant" |
      consequence_type == "downstream_gene_variant:missense_variant" |
      consequence_type == "exon_variant:missense_variant" |
      consequence_type == "missense_variant:upstream_gene_variant" |
      consequence_type == "missense_variant:intron_variant" |
      consequence_type == "upstream_gene_variant:missense_variant" |
      consequence_type == "missense_variant:3_prime_UTR_variant" |
      consequence_type == "missense_variant:downstream_gene_variant:exon_variant" |
      consequence_type == "missense_variant:exon_variant:downstream_gene_variant" |
      consequence_type == "3_prime_UTR_variant:missense_variant" |
      consequence_type == "downstream_gene_variant:missense_variant:upstream_gene_variant" |
      consequence_type == "downstream_gene_variant:missense_variant:exon_variant"
  )
synonymous_variant <-
  mutations_wCCF %>% filter(consequence_type == "synonymous_variant")

# define coding and non-coding mutations
non_coding <- rbind(intergenic_region, intergenic_region_mixed)
coding <- rbind(
  missense_variant,
  missense_variant_mixed,
  exon_variant,
  exon_variant_mixed,
  synonymous_variant
)

timing_kar <- function(ta) {
  ta2.0 <- ta[ta$karyotype == "2:0", ]
  ta2.0$timing <- ifelse(ta2.0$multiplicity == 1, "later", "early")
  ta2.1 <- ta[ta$karyotype == "2:1", ]
  ta2.1$timing <- ifelse(ta2.1$multiplicity == 1, "later", "early")
  ta2.2 <- ta[ta$karyotype == "2:2", ]
  ta2.2$timing <- ifelse(ta2.2$multiplicity == 1, "later", "early")
  
  ta3.0 <- ta[ta$karyotype == "3:0", ]
  ta3.0$timing <- ifelse(ta3.0$multiplicity <= 2, "later", "early")
  ta3.1 <- ta[ta$karyotype == "3:1", ]
  ta3.1$timing <- ifelse(ta3.1$multiplicity <= 2, "later", "early")
  ta3.2 <- ta[ta$karyotype == "3:2", ]
  ta3.2$timing <- ifelse(ta3.2$multiplicity <= 2, "later", "early")
  ta3.3 <- ta[ta$karyotype == "3:3", ]
  ta3.3$timing <- ifelse(ta3.3$multiplicity <= 2, "later", "early")
  
  ta4.0 <- ta[ta$karyotype == "4:0", ]
  ta4.0$timing <- ifelse(ta4.0$multiplicity <= 3, "later", "early")
  ta4.1 <- ta[ta$karyotype == "4:1", ]
  ta4.1$timing <- ifelse(ta4.1$multiplicity <= 3, "later", "early")
  ta4.2 <- ta[ta$karyotype == "4:2", ]
  ta4.2$timing <- ifelse(ta4.2$multiplicity <= 3, "later", "early")
  ta4.3 <- ta[ta$karyotype == "4:3", ]
  ta4.3$timing <- ifelse(ta4.3$multiplicity <= 3, "later", "early")
  
  ta5.0 <- ta[ta$karyotype == "5:0", ]
  ta5.0$timing <- ifelse(ta5.0$multiplicity <= 4, "later", "early")
  ta5.1 <- ta[ta$karyotype == "5:1", ]
  ta5.1$timing <- ifelse(ta5.1$multiplicity <= 4, "later", "early")
  ta5.2 <- ta[ta$karyotype == "5:2", ]
  ta5.2$timing <- ifelse(ta5.2$multiplicity <= 4, "later", "early")
  
  ta6.1 <- ta[ta$karyotype == "6:1", ]
  ta6.1$timing <- ifelse(ta6.1$multiplicity <= 5, "later", "early")
  ta6.2 <- ta[ta$karyotype == "6:2", ]
  ta6.2$timing <- ifelse(ta6.2$multiplicity <= 5, "later", "early")
  
  ta <-
    rbind(
      ta2.0,
      ta2.1,
      ta2.2,
      ta3.0,
      ta3.1,
      ta3.2,
      ta3.3,
      ta4.0,
      ta4.1,
      ta4.2,
      ta4.3,
      ta5.0,
      ta5.1,
      ta5.2
    )
  return(ta)
}

non_coding_timing <- timing_kar(non_coding)
coding_timing <- timing_kar(coding)

pie(table(non_coding_timing$timing) / length(non_coding_timing$timing))
pie(table(coding_timing$timing) / length(coding_timing$timing))

mutations_wCCF$multiplicity <-
  as.factor(mutations_wCCF$multiplicity)

# produce histogram of VAF over different karyotypes
pdf(file = paste0(
  "results/plots/000_paper_plots/00_Fig4",
  ifelse(clonal, "_clonal", "_noClonal"),
  ifelse(test, "_mockData", ""),
  ".pdf"
))
for (kar in c(
  "2:0",
  "2:1",
  "2:2",
  "3:0",
  "3:1",
  "3:2",
  "3:3",
  "4:0",
  "4:1",
  "4:2",
  "4:3",
  "5:0",
  "5:1",
  "5:2"
)) {
  t <- coding_timing %>% filter(karyotype == kar)
  
  p <- non_coding_timing %>% filter(karyotype == kar)
  
  par(mfrow = c(1, 2))
  print(
    ggplot(mutations_wCCF %>% filter(karyotype == kar)) +
      geom_histogram(
        aes(x = VAF / purity, fill = multiplicity),
        alpha = 0.3,
        binwidth = 0.01
      ) +
      theme_bw() +
      xlim(0, 1.3) +
      ggtitle(paste("PCAWG - karyotype", kar))
  )
  pie(table(t$timing) / sum(table(t$timing)),
      main = paste(kar, "- coding mutations"))
  pie(table(p$timing) / sum(table(p$timing)),
      main = paste(kar, "- non-coding mutations"))
  
}
dev.off()

rbind(table(non_coding_timing$timing),
      table(coding_timing$timing))[, 1:2]
fisher.test(rbind(
  table(non_coding_timing$timing),
  table(coding_timing$timing)
)[, 1:2], alternative = "greater")
chisq.test(rbind(
  table(non_coding_timing$timing),
  table(coding_timing$timing)
)[, 1:2])

source("https://www.r-statistics.com/wp-content/uploads/2010/02/Barnard.R.txt")

test_independence <- function(cod, ncod, kar) {
  cod <- cod %>% filter(karyotype == kar)
  ncod <- ncod %>% filter(karyotype == kar)
  
  t <- rbind(table(ncod$timing), table(cod$timing))[, 1:2]
  print(t)
  
  t_perc <- rbind(table(ncod$timing), table(cod$timing))[, 1:2]
  t_perc <-
    rbind(table(ncod$timing) / sum(table(ncod$timing)),
          table(cod$timing) / sum(table(cod$timing)))
  
  print(t_perc)
  n_mut <- dim(cod)[1] + dim(ncod)[1]
  return(append(fisher.test(t, alternative = "greater"), Barnard(t), n_mut))
}

ftot <- data.frame()
for (kar in c(
  "2:0",
  "2:1",
  "2:2",
  "3:0",
  "3:1",
  "3:2",
  "3:3",
  "4:0",
  "4:1",
  "4:2",
  "4:3",
  "5:0",
  "5:1",
  "5:2"
)) {
  print(kar)
  f <- test_independence(coding_timing, non_coding_timing, kar)
  
  ftot <-
    rbind(ftot,
          cbind(
            karyotype = kar,
            odds.ratio = f$estimate,
            p.value = f$p.value,
            n_muts = f[[8]]
          ))
}

ftot[, 2:4] <- apply(ftot[, 2:4], 2, as.numeric)
view(ftot)
write.table(file = "results/tables/04_timingAnalysis/statistics_coding_noncoding.txt")

timed_mutations <- rbind(coding_timing, non_coding_timing)
timed_mutations <-
  as.data.frame(table(timed_mutations$karyotypeF, timed_mutations$timing))

p <- ggplot(w) +
  geom_bar(aes(x = Var1, y = Freq, fill = Var2), stat = "summary") +
  theme_classic()

pdf(file = "results/plots/000_paper_plots/00_SupplementaryFig5.pdf")
print(p)
dev.off()

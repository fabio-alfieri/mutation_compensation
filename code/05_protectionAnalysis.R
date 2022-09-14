# GENE-PROTECTION ANALYSIS
# This script computes protection index for each gene and thus perform enrichment
# analysis (Gene Ontology) on both protected and unprotected gene sets

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!require("clusterProfiler", quietly = TRUE)){
  library(BiocManager)
  BiocManager::install("clusterProfiler")
}
suppressMessages({
  library(clusterProfiler)
  library(readxl)
  library(tibble)
  library(ggplot2)
  library(tidyr)
  library(readr)
  library(stringr)
  library(data.table)
  library(dplyr)
  library(ggExtra)
  library(purrr)
  library(clusterProfiler)
  library("org.Hs.eg.db")
  library(AnnotationDbi)
})

setwd("../")

# selected tumor types that showed a positive and significant correlation estimates
# based on Spearman's (cor: amplifications ~ all_mutations)
tumor_types <- c(  
  "LUAD", "LUSC", "BRCA", "CESC", "THCA", "HNSC", "PAAD", "GBMLGG", "COADREAD"
)

cat("\n\n >> Calculate protection index \n\n")
statistics <- data.frame()
for (tumor_type in tumor_types) {
  cat("\t\t", tumor_type, "\n")
  
  suppressMessages({
  genelevel <-
    read.csv(paste0(
      "results/tables/00_geneLevelAnalysis/",
      tumor_type,
      "_geneLevel.csv"
    ))
  })
  
  # add a psuedocount for mutations
  genelevel$MutationsDiploid_norm <-
    genelevel$MutationsDiploid_norm + quantile(genelevel$MutationsDiploid_norm,
                                               prob = c(0:100 / 100))[quantile(genelevel$MutationsDiploid_norm,
                                                                               prob = c(0:100 / 100)) > 0][1]
  
  genelevel <-
    genelevel[!duplicated(paste0(genelevel$Gene, genelevel$GeneLegth)),]
  genelevel <-
    genelevel[!(
      is.na(genelevel$Amplification) |
        is.na(genelevel$MutationsDiploid_norm) |
        is.na(genelevel$GeneTPM)
    ),]
  
  # remove very lowly expressed/or no expressed genes
  genelevel <-
    genelevel[genelevel$GeneTPM >= quantile(genelevel$GeneTPM, prob = c(0:100 /
                                                                          100)[5]), ]
  
  # calcualte the mutation score
  genelevel$mutation.score <-
    1 - log10(genelevel$MutationsDiploid_norm) / min(log10(genelevel$MutationsDiploid_norm))
  
  # define the protection index (PI) as amplifications - mutations
  genelevel$PI <-
    genelevel$Amplification - (genelevel$mutation.score * 2)
  
  genelevel$PI_norm <- genelevel$PI/(sum(genelevel$Amplification)+2*sum(genelevel$mutation.score))
  
  # define protected and unprotected gene sets according to PI
  protected <- data.frame()
  unprotected <- data.frame()
  for (chr in 1:22) {
    protected <-
      rbind(protected, cbind(genelevel[genelevel$Chromosome == chr &
                               genelevel$PI >= quantile(genelevel[genelevel$Chromosome == chr, ]$PI, prob = 0.9), ]$Gene,
                             genelevel[genelevel$Chromosome == chr &
                                         genelevel$PI >= quantile(genelevel[genelevel$Chromosome == chr, ]$PI, prob = 0.9), ]$PI,
                             genelevel[genelevel$Chromosome == chr &
                                         genelevel$PI >= quantile(genelevel[genelevel$Chromosome == chr, ]$PI, prob = 0.9), ]$PI_norm))
    unprotected <-
      rbind(unprotected, cbind(genelevel[genelevel$Chromosome == chr &
                                 genelevel$PI <= quantile(genelevel[genelevel$Chromosome == chr, ]$PI, prob = 0.05), ]$Gene,
                               genelevel[genelevel$Chromosome == chr &
                                           genelevel$PI <= quantile(genelevel[genelevel$Chromosome == chr, ]$PI, prob = 0.05), ]$PI,
                               genelevel[genelevel$Chromosome == chr &
                                           genelevel$PI <= quantile(genelevel[genelevel$Chromosome == chr, ]$PI, prob = 0.05), ]$PI_norm))
  }
  
  genelevel$protection.status <-
    ifelse(
      genelevel$Gene %in% protected$V1,
      "PROTECTED",
      ifelse(genelevel$Gene %in% unprotected$V2, "UNPROTECTED", "ns")
    )
  
  # define the background gene set as the union of protected and unprotected gene sets
  
  # output of the analysis - write tables
  system("mkdir -p results/tables/00_geneLevelAnalysis/tmp/")
  
  protected[,2] <- as.numeric(protected[,2])
  protected <- protected[order(protected$V2, decreasing = T),]
  unprotected[,2] <- as.numeric(unprotected[,2])
  unprotected <- unprotected[order(unprotected$V2, decreasing = T),]
  
  
  write_tsv(
    protected,
    # quote = F,
    # row.names = F,
    file = paste0(
      "results/tables/00_geneLevelAnalysis/tmp/",
      tumor_type,
      "_protected.tsv"
    )
  )
  write_tsv(
    unprotected,
    # quote = F,
    # row.names = F,
    file = paste0(
      "results/tables/00_geneLevelAnalysis/tmp/",
      tumor_type,
      "_unprotected.tsv"
    )
  )
}

# define a common list of protected and unprotected according the overlap between
# analyzed cancer types

protected <- data.frame()
unprotected <- data.frame()
background <- data.frame()

protected_all <- data.frame()
unprotected_all <- data.frame()
background_all <- data.frame()

for (tumor_type in tumor_types) {
  suppressMessages({
  protected <-
    read_tsv(
      paste0(
        "results/tables/00_geneLevelAnalysis/tmp/",
        tumor_type,
        "_protected.tsv"
      )
    )
  unprotected <-
    read_tsv(
      paste0(
        "results/tables/00_geneLevelAnalysis/tmp/",
        tumor_type,
        "_unprotected.tsv"
      )
    )
  })
  
  background <- rbind(protected, unprotected)
  
  protected$tumor_type <- tumor_type
  unprotected$tumor_type <- tumor_type
  background$tumor_type <- tumor_type
  
  protected_all <- rbind(protected_all, protected)
  unprotected_all <- rbind(unprotected_all, unprotected)
  background_all <- rbind(background_all, background)
}

rm(tumor_type)
rm(tumor_types)
rm(protected)
rm(background)
rm(unprotected)

cat("\n\n >> Extacting common protected and unprotected genes across tumor types \n\n")
calculate_mean <- function(df, n){
  df1 <- df[df$V1 %in% names(table(df[,1])[table(df[,1]) >= n]),]
  df2 <- data.frame()
  for(gene in levels(factor(df1$V1))){
    df2 <- rbind(df2, cbind(gene, mean = mean(df1[df1[,1] == gene,]$V3)))
  }
  return(df2)
}

write_tsv(
  calculate_mean(protected_all, 3),
  "results/tables/00_geneLevelAnalysis/00_protected.tsv"
  # # quote = F,
  # row.names = F
)
write_tsv(
  calculate_mean(unprotected_all, 3),
  "results/tables/00_geneLevelAnalysis/00_unprotected.tsv"
  # quote = F,
  # row.names = F
)
write_tsv(
  calculate_mean(background_all, 3),
  "results/tables/00_geneLevelAnalysis/00_background.tsv"
  # quote = F,
  # row.names = F
)

rm(genelevel)
rm(background_all)
rm(protected_all)
rm(unprotected_all)

# from CCLE (Cancer Cell Line Enciclopedia) https://depmap.org/portal/download/
cat("\n\n >> Calculating essentiality scores for each gene of the two categories \n\n")
load(file = "data/misc/mean_crispr_effect.RData")
protected <-
  read.table("results/tables/00_geneLevelAnalysis/00_protected.tsv")
unprotected <-
  read.table("results/tables/00_geneLevelAnalysis/00_unprotected.tsv")

t.test(mean_crispr_effect[mean_crispr_effect$genes %in% protected$V1,]$mean,
       mean_crispr_effect[mean_crispr_effect$genes %in% unprotected$V1,]$mean,
       alternative = "less")
w <-
  wilcox.test(mean_crispr_effect[mean_crispr_effect$genes %in% protected$V1,]$mean,
              mean_crispr_effect[mean_crispr_effect$genes %in% unprotected$V1,]$mean,
              alternative = "less")

cat(" \n > Producing Fig. 5 with and without outliers \n\n")
pdf("results/plots/000_paper_plots/00_Fig5c.pdf")
boxplot(mean_crispr_effect[mean_crispr_effect$genes %in% protected$V1,]$mean,
        mean_crispr_effect[mean_crispr_effect$genes %in% unprotected$V1,]$mean)
boxplot(mean_crispr_effect[mean_crispr_effect$genes %in% protected$V1,]$mean,
        mean_crispr_effect[mean_crispr_effect$genes %in% unprotected$V1,]$mean,
        outline = F)
abline(h = 0, lty = "dotted")
dev.off()

print(w)


# GO analysis on protected gene set ----
cat("\n\n >> Gene Ontology analysis \n\n")
cat(" > Analyzing the PROTECTED gene set \n\n")
library(enrichplot)
entrezid_buffered <-
  mapIds(org.Hs.eg.db, protected$V1, 'ENTREZID', 'SYMBOL')
entrezid_unbuffered <-
  mapIds(org.Hs.eg.db, unprotected$V1, 'ENTREZID', 'SYMBOL')
entrezid_background <-
  mapIds(org.Hs.eg.db,
         c(protected$V1, unprotected$V1),
         'ENTREZID',
         'SYMBOL')

eGO_bufferedBP <-
  clusterProfiler::enrichGO(
    entrezid_buffered,
    'org.Hs.eg.db',
    ont = "ALL",
    universe = entrezid_background,
    readable = TRUE,
    pAdjustMethod = "fdr",
    minGSSize = 2
  )
head(eGO_bufferedBP, n = 100)[c(1,3,6,7,10)]

# UNNPROTECTED -----
cat("\n\n > Analyzing the PROTECTED gene set \n\n")
eGO_unbufferedBP <-
  clusterProfiler::enrichGO(
    entrezid_unbuffered,
    'org.Hs.eg.db',
    ont = "ALL",
    universe = entrezid_background,
    readable = TRUE,
    pAdjustMethod =  "fdr",
    minGSSize = 2
  )
head(eGO_unbufferedBP, n = 100)[c(1,3,6,7,10)]


# filter results based on qvalue and GeneCount
buffered_results <- eGO_bufferedBP@result
buffered_results <-
  buffered_results[buffered_results$p.adjust <= 0.05, ]

unbuffered_results <- eGO_unbufferedBP@result

unbuffered_results <-
  unbuffered_results[unbuffered_results$p.adjust <= 0.05, ]

write_tsv(buffered_results,
          file = "results/tables/00_geneLevelAnalysis/00_GOprotected.tsv")
write_tsv(unbuffered_results,
          file = "results/tables/00_geneLevelAnalysis/00_GOunprotected.tsv")

write.table(
  noquote(cbind(
    buffered_results$ID, buffered_results$p.adjust
  )),
  file = "results/tables/00_geneLevelAnalysis/01_protected_4REVIGO.txt",
  quote = F,
  row.names = F,
  col.names = F
)
write.table(
  noquote(cbind(
    unbuffered_results$ID, unbuffered_results$p.adjust
  )),
  file = "results/tables/00_geneLevelAnalysis/01_unprotected_4REVIGO.txt",
  quote = F,
  row.names = F,
  col.names = F
)

# now go on http://revigo.irb.hr/ to plot results

# comparison with CRISPR common essential and non-essential genes
suppressMessages({
  common_essential <-
    read_csv(file = "data/misc/CRISPR_common_essentials.csv")
  nonessential <-
    read_csv(file = "data/misc/CRISPR_nonessentials_NegativeCTRs.csv")
})
table(protected$V1 %in% common_essential$gene)
table(protected$V1 %in% nonessential$gene)

cat("\n\n >> Contingency tables \n\n")
rbind(protected = cbind(table(protected$V1 %in% common_essential$gene)[2],
                        table(protected$V1 %in% nonessential$gene)[2]),
      unprotected = cbind(table(unprotected$V1 %in% common_essential$gene)[2],
                          table(unprotected$V1 %in% nonessential$gene)[2]))
f <- fisher.test(rbind(protected = cbind(table(protected$V1 %in% common_essential$gene)[2],
                        table(protected$V1 %in% nonessential$gene)[2]),
      unprotected = cbind(table(unprotected$V1 %in% common_essential$gene)[2],
                          table(unprotected$V1 %in% nonessential$gene)[2])))

rbind(
  table(protected$V1 %in% nonessential$gene),
  table(unprotected$V1 %in% nonessential$gene)
)
fisher.test(rbind(
  table(protected$V1 %in% nonessential$gene),
  table(unprotected$V1 %in% nonessential$gene)
))

 rbind(
  table(protected$V1 %in% common_essential$gene),
  table(unprotected$V1 %in% common_essential$gene)
)
fisher.test(rbind(
  table(protected$V1 %in% common_essential$gene),
  table(unprotected$V1 %in% common_essential$gene)
))


cat("\n\n OUTPUT of the script: \n \t (1) raw tables path: results/tables/00_geneLevelAnalysis/ \n")
cat("\t (2) paper-like plots path: results/plots/000_paper_plots/ \n\n")

# GENE-PROTECTION ANALYSIS
#
# This script computes protection index for each gene and thus perform enrichment
# analysis (Gene Ontology) on both protected and unprotected gene sets

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!require("clusterProfiler", quietly = TRUE)){
  library(BiocManager)
  BiocManager::install("clusterProfiler")
}
if (!require("org.Hs.eg.db", quietly = TRUE)){
  BiocManager::install("org.Hs.eg.db")
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
  library(ggpubr)
  library(purrr)
  library(clusterProfiler)
  library("org.Hs.eg.db")
  library(AnnotationDbi)
  library(parallel)
  library(optparse)
})

option_list = list(
  make_option(
    c("-p", "--runPermutations"),
    type = "character",
    default = "n",
    help = "Options are: ([y]/n)
                It requires some minutes (up to 1 h). 
                It needs parallelization",
    metavar = ""
  )
)

opt_parser = OptionParser(option_list = option_list)

opt = parse_args(opt_parser)


if (opt$runPermutations == "n") {
  cat("\n\n >> You chose default options: \n\t --runPermutations = 'n' \n\n")
}

if (!any(opt$runPermutations %in% c("y", "n"))) {
  print_help(opt_parser)
  stop("typo in the analysis flag, plase see above for the available options!",
       call. = FALSE)
}


setwd("../")
# setwd("~/Desktop_linux/mutation_compensation")

OGs <-read.csv(file = "data/CancerGenes/OG_list.tsv", sep = "\t")
TSGs <-read.csv(file = "data/CancerGenes/TSG_list.tab", sep = "\t")

# selected tumor types that showed a positive and significant correlation estimates
# based on Spearman's (cor: amplifications ~ all_mutations)
tumor_types <- c(
  "LUAD",
  "LUSC",
  "BRCA",
  "CESC",
  "THCA",
  "HNSC",
  "PAAD",
  "COADREAD",
  "GBMLGG"
)


# produce permutation on Pi scores
if(opt$runPermutations == "y"){
  cat("\n\n >> Run permutations on Protection index \n\n")
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
    
    # calculate the mutation score
    genelevel$mutation.score <-
      1 - log10(genelevel$MutationsDiploid_norm) / min(log10(genelevel$MutationsDiploid_norm))
    
    # define the protection index (PI) as amplifications - mutations
    genelevel$PI <-
      genelevel$Amplification - (2*genelevel$mutation.score)
    
    # genelevel$PI_norm <- genelevel$PI/(sum(genelevel$Amplification)+2*sum(genelevel$mutation.score))
    
    for(i in 1:10000){
      set.seed(i)
      genelevel[,ncol(genelevel)+1] <- sample(genelevel$Amplification) - sample(2*genelevel$mutation.score)
    }
    
    p_distributions <- mclapply(1:nrow(genelevel), mc.cores = 35, function(x){
      vect_cum_dist <- cume_dist(as.numeric(unlist(c(genelevel[x,19:ncol(genelevel)]))))
      vect_cum_dist <- as.data.frame(vect_cum_dist)
      vect_cum_dist$PI_scores <- as.numeric(unlist(c(genelevel[x,19:ncol(genelevel)])))
      res = vect_cum_dist[vect_cum_dist$PI_scores == genelevel[x,"PI"],][1,1]
      names(res) <- genelevel[x,"Gene"]
      return(res)
    })
    
    p_distributions <- as.data.frame(unlist(p_distributions))
    p_distributions$Gene <- rownames(p_distributions)
    rownames(p_distributions) <- NULL
    p_distributions$protected.pvals <- 1 - p_distributions$`unlist(p_distributions)`
    colnames(p_distributions) <- c("unprotected.pdistr", "Gene", "protected.pdistr")
    p_distributions <- p_distributions[,c(2,1,3)]
    
    genes[[tumor_type]] <- p_distributions
  }
  
  saveRDS(genes, file = "data/misc/Pi_permutations.rds")
}


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
  
  # colnames(genelevel)[10] <- "MutationsDiploid_norm"
  
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
  
  # calculate the mutation score
  genelevel$mutation.score <-
    1 - log10(genelevel$MutationsDiploid_norm) / min(log10(genelevel$MutationsDiploid_norm))
  
  # define the protection index (PI) as amplifications - mutations
  genelevel$PI <-
    genelevel$Amplification - (genelevel$mutation.score * 2)
  
  genelevel$PI_norm <- genelevel$PI/(sum(genelevel$Amplification)+2*sum(genelevel$mutation.score))
  
  genes_pdist <- readRDS("data/misc/Pi_permutations.rds")
  genes_pdist_tissue <- genes_pdist[[tumor_type]]
  
  genelevel <- left_join(genes_pdist_tissue, genelevel, by = "Gene")
  
  genelevel <- genelevel[!genelevel$Gene %in% unlist(strsplit(OGs$Gene.names, " ")),]
  genelevel <- genelevel[!genelevel$Gene %in% unlist(strsplit(TSGs$Gene.names, " ")),]
  
  # define protected and unprotected gene sets according to PI
  protected <- data.frame()
  unprotected <- data.frame()
  for (chr in 1:22) {
    protected <-
      rbind(protected, cbind(genelevel[genelevel$Chromosome == chr &
                               genelevel$PI >= quantile(genelevel[genelevel$Chromosome == chr, ]$PI, prob = 0.94), ]$Gene,
                             genelevel[genelevel$Chromosome == chr &
                                         genelevel$PI >= quantile(genelevel[genelevel$Chromosome == chr, ]$PI, prob = 0.94), ]$PI,
                             genelevel[genelevel$Chromosome == chr &
                                         genelevel$PI >= quantile(genelevel[genelevel$Chromosome == chr, ]$PI, prob = 0.94), ]$PI_norm,
                             genelevel[genelevel$Chromosome == chr &
                                         genelevel$PI >= quantile(genelevel[genelevel$Chromosome == chr, ]$PI, prob = 0.94), ]$protected.pdistr))
    unprotected <-
      rbind(unprotected, cbind(genelevel[genelevel$Chromosome == chr &
                                 genelevel$PI <= quantile(genelevel[genelevel$Chromosome == chr, ]$PI, prob = 0.06), ]$Gene,
                               genelevel[genelevel$Chromosome == chr &
                                           genelevel$PI <= quantile(genelevel[genelevel$Chromosome == chr, ]$PI, prob = 0.06), ]$PI,
                               genelevel[genelevel$Chromosome == chr &
                                           genelevel$PI <= quantile(genelevel[genelevel$Chromosome == chr, ]$PI, prob = 0.06), ]$PI_norm,
                               genelevel[genelevel$Chromosome == chr &
                                           genelevel$PI <= quantile(genelevel[genelevel$Chromosome == chr, ]$PI, prob = 0.06), ]$unprotected.pdistr))
  }
  
  genelevel$protection.status <-
    ifelse(
      genelevel$Gene %in% protected$V1,
      "PROTECTED",
      ifelse(genelevel$Gene %in% unprotected$V1, "UNPROTECTED", "ns")
    )
  
  # define the background gene set as the union of protected and unprotected gene sets
  
  # output of the analysis - write tables
  system("mkdir -p results/tables/00_geneLevelAnalysis/tmp/")
  
  protected[,2:4] <- apply(protected[,2:4], 2, as.numeric)
  protected <- protected[order(protected$V2, decreasing = T),]
  unprotected[,2:4] <- apply(unprotected[,2:4], 2, as.numeric)
  unprotected <- unprotected[order(unprotected$V2, decreasing = T),]
  
  colnames(protected) <- c("Gene", "PI", "PI_norm", "p.distr")
  colnames(unprotected) <- c("Gene", "PI", "PI_norm", "p.distr")
  
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
  df1 <- df[df$Gene %in% names(table(df[,1])[table(df[,1]) >= n]),]
  df2 <- data.frame()
  for(gene in levels(factor(df1$Gene))){
    df2 <- rbind(df2, cbind(gene, 
                            mean_PInorm = mean(df1[df1[,1] == gene,]$PI_norm), 
                            mean_PI = mean(df1[df1[,1] == gene,]$PI),
                            p.distr = mean(df1[df1[,1] == gene,]$p.distr)))
  }
  return(df2)
}

common <- 3
write_tsv(
  calculate_mean(protected_all, common),
  "results/tables/00_geneLevelAnalysis/00_protected.tsv"
)
write_tsv(
  calculate_mean(unprotected_all, common),
  "results/tables/00_geneLevelAnalysis/00_unprotected.tsv"
)

protected <-
  read.table("results/tables/00_geneLevelAnalysis/00_protected.tsv", header = T)
unprotected <-
  read.table("results/tables/00_geneLevelAnalysis/00_unprotected.tsv", header = T)

protected <- protected %>% filter(p.distr <= 0.1)
unprotected <- unprotected %>% filter(p.distr <= 0.1)

rm(genelevel)
rm(background_all)
rm(protected_all)
rm(unprotected_all)

# from CCLE (Cancer Cell Line Enciclopedia) https://depmap.org/portal/download/
cat("\n\n >> Calculating essentiality scores for each gene of the two categories \n\n")
load(file = "data/misc/mean_crispr_effect.RData")

t.test(mean_crispr_effect[mean_crispr_effect$genes %in% protected$gene,]$mean,
       mean_crispr_effect[mean_crispr_effect$genes %in% unprotected$gene,]$mean,
       alternative = "less")
(w <-
  wilcox.test(mean_crispr_effect[mean_crispr_effect$genes %in% protected$gene,]$mean,
              mean_crispr_effect[mean_crispr_effect$genes %in% unprotected$gene,]$mean))

cat(" \n > Producing Fig. 5 with and without outliers \n\n")
pdf("results/plots/000_paper_plots/00_Fig5c.pdf")
boxplot(mean_crispr_effect[mean_crispr_effect$genes %in% protected$gene,]$mean,
        mean_crispr_effect[mean_crispr_effect$genes %in% unprotected$gene,]$mean)
boxplot(mean_crispr_effect[mean_crispr_effect$genes %in% protected$gene,]$mean,
        mean_crispr_effect[mean_crispr_effect$genes %in% unprotected$gene,]$mean,
        outline = F)
abline(h = 0, lty = "dotted")
dev.off()

# GO analysis on protected gene set ----
cat("\n\n >> Gene Ontology analysis \n\n")
cat(" > Analyzing the PROTECTED gene set \n\n")
library(enrichplot)
entrezid_buffered <-
  mapIds(org.Hs.eg.db, protected$gene, 'ENTREZID', 'SYMBOL')
entrezid_unbuffered <-
  mapIds(org.Hs.eg.db, unprotected$gene, 'ENTREZID', 'SYMBOL')
entrezid_background <-
  mapIds(org.Hs.eg.db,
         c(protected$gene, unprotected$gene),
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

# buffered rrvgo
library(rrvgo)
revigo_buff <- noquote(cbind(
  eGO_bufferedBP$ID, eGO_bufferedBP$p.adjust))
GO_terms <- data.frame()
for(cat in c("BP", "MF", "CC")){
  simMatrix <- calculateSimMatrix(revigo_buff[,1],
                                  orgdb = "org.Hs.eg.db",
                                  ont = cat,
                                  method = "Rel")
  scores <- setNames(-log10(as.numeric(revigo_buff[,2])), revigo_buff[,1])
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
  GO_terms_tmp <- reducedTerms
  GO_terms_tmp$cat <- cat
  GO_terms <- rbind(GO_terms, GO_terms_tmp)
}
treemapPlot(GO_terms)

pdf(file = "results/plots/00_geneLevelAnalysis/buffered_DNARNA_BP.pdf",
    width = 15, height = 15) 
filter(GO_terms, cat == "BP") %>% treemapPlot()
dev.off()
pdf(file = "results/plots/00_geneLevelAnalysis/buffered_DNARNA_MF.pdf",
    width = 15, height = 15) 
filter(GO_terms, cat == "MF") %>% treemapPlot()
dev.off()
pdf(file = "results/plots/00_geneLevelAnalysis/buffered_DNARNA_CC.pdf", 
    width = 15, height = 15) 
filter(GO_terms, cat == "CC") %>% treemapPlot()
dev.off()


# UNNPROTECTED -----
cat("\n\n > Analyzing the UNPROTECTED gene set \n\n")
eGO_unbufferedBP <-
  clusterProfiler::enrichGO(
    entrezid_unbuffered,
    'org.Hs.eg.db',
    ont = "ALL",
    readable = TRUE,
    pAdjustMethod =  "fdr",
    minGSSize = 20
  )
head(eGO_unbufferedBP, n = 100)[c(1,3,6,7,10)]

# unbuffered rrvgo
library(rrvgo)
revigo_unbuff <- noquote(cbind(
  eGO_unbufferedBP$ID, eGO_unbufferedBP$p.adjust))
GO_terms <- data.frame()
for(cat in c("BP", "MF", "CC")){
  simMatrix <- calculateSimMatrix(revigo_unbuff[,1],
                                  orgdb = "org.Hs.eg.db",
                                  ont = cat,
                                  method = "Rel")
  scores <- setNames(-log10(as.numeric(revigo_unbuff[,2])), revigo_unbuff[,1])
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
  GO_terms_tmp <- reducedTerms
  GO_terms_tmp$cat <- cat
  GO_terms <- rbind(GO_terms, GO_terms_tmp)
}
treemapPlot(GO_terms)

pdf(file = "results/plots/00_geneLevelAnalysis/unbuffered_DNARNA_BP.pdf",
    width = 15, height = 15) 
filter(GO_terms, cat == "BP") %>% treemapPlot()
dev.off()
pdf(file = "results/plots/00_geneLevelAnalysis/unbuffered_DNARNA_MF.pdf",
    width = 15, height = 15) 
filter(GO_terms, cat == "MF") %>% treemapPlot()
dev.off()
pdf(file = "results/plots/00_geneLevelAnalysis/unbuffered_DNARNA_CC.pdf", 
    width = 15, height = 15) 
filter(GO_terms, cat == "CC") %>% treemapPlot()
dev.off()


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

library("openxlsx")
write.xlsx(buffered_results,
                file = "results/tables/00_geneLevelAnalysis/00_GOprotected.xlsx")
write.xlsx(unbuffered_results,
                file = "results/tables/00_geneLevelAnalysis/00_GOunprotected.xlsx")

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
table(protected$gene %in% common_essential$gene)
table(protected$gene %in% nonessential$gene)

cat("\n\n >> Contingency tables \n\n")
rbind(protected = cbind(table(protected$gene %in% common_essential$gene)[2],
                        table(protected$gene %in% nonessential$gene)[2]),
      unprotected = cbind(table(unprotected$gene %in% common_essential$gene)[2],
                          table(unprotected$gene %in% nonessential$gene)[2]))

f <- fisher.test(rbind(protected = cbind(table(protected$gene %in% common_essential$gene)[2],
                        table(protected$gene %in% nonessential$gene)[2]),
      unprotected = cbind(table(unprotected$gene %in% common_essential$gene)[2],
                          table(unprotected$gene %in% nonessential$gene)[2])))

rbind(
  table(protected$gene %in% nonessential$gene),
  table(unprotected$gene %in% nonessential$gene)
)
fisher.test(rbind(
  table(protected$gene %in% nonessential$gene),
  table(unprotected$gene %in% nonessential$gene)
))
 
rbind(
  table(protected$gene %in% common_essential$gene),
  table(unprotected$gene %in% common_essential$gene)
)
fisher.test(rbind(
  table(protected$gene %in% common_essential$gene),
  table(unprotected$gene %in% common_essential$gene)
))


protected_ess <- protected[protected$gene %in% common_essential$gene,]$gene
protected_noness <- protected[protected$gene %in% nonessential$gene,]$gene
unprotected_ess <- unprotected[unprotected$gene %in% common_essential$gene,]$gene
unprotected_noness <- unprotected[unprotected$gene %in% nonessential$gene,]$gene

write.xlsx(as.data.frame(protected_ess),
      file = "results/tables/00_geneLevelAnalysis/00_protected_ess.xlsx")
write.xlsx(as.data.frame(protected_noness), 
           file = "results/tables/00_geneLevelAnalysis/00_protected_noness.xlsx")
write.xlsx(as.data.frame(unprotected_ess),
           file = "results/tables/00_geneLevelAnalysis/00_unprotected_ess.xlsx")
write.xlsx(as.data.frame(unprotected_noness),
           file = "results/tables/00_geneLevelAnalysis/00_unprotected_noness.xlsx")

## check expression of protected and unprotected genes ----
pdf(file = "results/plots/00_geneLevelAnalysis/genexpression.pdf")
for(tumor_type in tumor_types){
  suppressMessages({
    genelevel <-
      read.csv(paste0(
        "results/tables/00_geneLevelAnalysis/",
        tumor_type,
        "_geneLevel.csv"
      ))
  })
  protected_tpm <- genelevel %>% filter(Gene %in% protected$gene)
  unprotected_tpm <- genelevel %>% filter(Gene %in% unprotected$gene)
  tpm <- as.data.frame(rbind(cbind(log10(as.numeric(protected_tpm$GeneTPM)), "protected"),
                             cbind(log10(as.numeric(unprotected_tpm$GeneTPM)), "unprotected")))
  tpm[,1] <- as.numeric(tpm[,1])
  tpm[,2] <- as.factor(tpm[,2])
  
  png(file = paste0("results/plots/00_geneLevelAnalysis/",tumor_type,"_genexpression.png"))
  print(ggboxplot(tpm, x = "V2", y = "V1",
            palette = "jco") + 
    stat_compare_means() +
    ggtitle(tumor_type))
  dev.off()
}
dev.off()

cat("\n\n OUTPUT of the script: \n \t (1) raw tables path: results/tables/00_geneLevelAnalysis/ \n")
cat("\t (2) paper-like plots path: results/plots/000_paper_plots/ \n\n")

# BRCA-PAAD analysis
#
# Compensatory events partially explain tissue-specific amplification patterns

cat("\n\n >> Loading libraries \n\n")
suppressMessages({
  require(ggplot2)
  require(stringr)
  library(readr)
  require(mgcv)
  require(ggpubr)
  library(readxl)
  library(optparse)
  library(dplyr)
  library(utils)
  library(tidyverse)
})

setwd("../")

averageChrCoverage <-
  read.table("data/FireBrowse_CNAs/averageChrCoverage.txt")
averageChrCoverage <- averageChrCoverage[c(1:22),]

length_bin_CNA <- 1000000

cat(" \n > Calculate correlations between mu score and amplification frequency at 36Mbp\n\n")
genes_in_bin <- data.frame()
suppressMessages({
  for (chr in 1:22) {
    bin_gene <-
      read.table(
        paste0(
          "data/ChromosomeGeneStructure/chr_",
          chr,
          "_binSize_",
          length_bin_CNA,
          ".txt"
        )
      )
    bin_gene <- bin_gene[bin_gene$gene_count != 0,]
    
    bin_gene$bin_merged <-
      rep(1:round(dim(bin_gene)[1] / 36 + 0.5), each = round(dim(bin_gene)[1] /
                                                               round(dim(bin_gene)[1] / 36 + 0.5) + 0.5))[1:nrow(bin_gene)]
    
    bin_gene <-
      bin_gene[!bin_gene$bin_end <= averageChrCoverage[averageChrCoverage$chr == chr,]$start,]
    bin_gene <-
      bin_gene[!bin_gene$bin_end >= averageChrCoverage[averageChrCoverage$chr == chr,]$end,]
    
    for (i in 1:dim(bin_gene)[1]) {
      genes_in_bin <- rbind(
        genes_in_bin,
        cbind(
          chr = chr,
          bin = bin_gene[i,]$bin[1],
          gene = unlist(strsplit(bin_gene[i,]$gene_id, "/")),
          bin_id = paste0(chr, "_", bin_gene[i,]$bin),
          bin_merged36 = bin_gene[i,]$bin_merged,
          bin_id_merged36 = paste0(chr, "_", bin_gene[i,]$bin_merged)
        )
      )
    }
  }
})

# write.table(as.data.frame(genes_in_bin), file = "./bin_gene_position.tsv", sep = "\t")
rm(averageChrCoverage)
rm(bin_gene)

# flag both FALSE if you need the main analysis
CTR <- TRUE
if (CTR) {
  stat_estimates_CTR <- data.frame()
  stat_signif_CTR <- data.frame()
}

tumor_types <- c("BRCA", "PAAD")

stat_estimates <- data.frame()
stat_signif <- data.frame()

cat(" \n > Calculate OG and GO scores and produce Fig. S3\n\n")
suppressWarnings({
  suppressMessages({
    for (tumor_type in tumor_types) {
      merged_1 <- 36
      bin_values <-
        read.table(
          file = paste0(
            "results/tables/02_produceStatistics/amplifications_",
            tumor_type,
            "_",
            merged_1,
            "Mbp_table.txt"
          )
        )
      
      bin_values$bin_id_merged36 <-
        paste0(bin_values$chr, "_", bin_values$resize)
      
      # https://www.sciencedirect.com/science/ article/pii/S0092867413012877?via%3Dihub
      # Table 4 - TUSON Explorer Prediction of TSGs and OGs on Single Tumor Types
      tuson_OGs <-
        read_xlsx("data/misc/Davoli2013_TUSON.xlsx", sheet = 2)
      tuson_TSGs <-
        read_xlsx("data/misc/Davoli2013_TUSON.xlsx", sheet = 1)
      
      
      nOGs <- 250
      nTSGs <- 250
      genes_in_bin$OG.status.pancancer <-
        ifelse(genes_in_bin$gene %in% tuson_OGs[order(tuson_OGs$`PAN-Cancer_p-value`, decreasing = F),]$Gene[1:nOGs], "OG", "noOG")
      genes_in_bin$TSG.status.pancancer <-
        ifelse(genes_in_bin$gene %in% tuson_TSGs[order(tuson_TSGs$`PAN-Cancer_p-value`, decreasing = F),]$Gene[1:nTSGs], "TSG", "noTSG")
      
      genes_in_bin$value.OGs <- NULL
      genes_in_bin$value.TSGs <- NULL
      
      
      if (tumor_type == "BRCA") {
        genes_in_bin <- genes_in_bin[!is.na(genes_in_bin$chr),]
        
        genes_in_bin$OG.status <-
          ifelse(genes_in_bin$gene %in% tuson_OGs[order(tuson_OGs$`Breast_p-value`, decreasing = F),]$Gene[1:nOGs], "OG", "noOG")
        genes_in_bin$TSG.status <-
          ifelse(genes_in_bin$gene %in% tuson_TSGs[order(tuson_TSGs$`Breast_p-value`, decreasing = F),]$Gene[1:nTSGs],
                 "TSG",
                 "noTSG")
        
        tuson_OGs <-
          tuson_OGs[order(tuson_OGs$`Breast_p-value`, decreasing = F),]
        tuson_TSGs <-
          tuson_TSGs[order(tuson_TSGs$`Breast_p-value`, decreasing = F),]
        
        tuson_OGs$value.OGs <-
          c(nOGs:1, rep(0, length(tuson_OGs$Gene) - nOGs))
        tuson_TSGs$value.TSGs <-
          c(nTSGs:1, rep(0, length(tuson_TSGs$Gene) - nTSGs))
        
        colnames(tuson_OGs)[1] <- "gene"
        colnames(tuson_TSGs)[1] <- "gene"
        
        genes_in_bin <-
          left_join(genes_in_bin, tuson_OGs[, c(1, 44)], by = "gene")
        genes_in_bin <-
          left_join(genes_in_bin, tuson_TSGs[, c(1, 44)], by = "gene")
        
      } else if (tumor_type == "PAAD") {
        genes_in_bin <- genes_in_bin[!is.na(genes_in_bin$chr),]
        
        genes_in_bin$OG.status <-
          ifelse(genes_in_bin$gene %in% tuson_OGs[order(tuson_OGs$`Pancreas_p-value`, decreasing = F),]$Gene[1:nOGs], "OG", "noOG")
        genes_in_bin$TSG.status <-
          ifelse(genes_in_bin$gene %in% tuson_TSGs[order(tuson_TSGs$`Pancreas_p-value`, decreasing = F),]$Gene[1:nTSGs],
                 "TSG",
                 "noTSG")
        
        tuson_OGs <-
          tuson_OGs[order(tuson_OGs$`Pancreas_p-value`, decreasing = F),]
        tuson_TSGs <-
          tuson_TSGs[order(tuson_TSGs$`Pancreas_p-value`, decreasing = F),]
        
        tuson_OGs$value.OGs <-
          c(nOGs:1, rep(0, length(tuson_OGs$Gene) - nOGs))
        tuson_TSGs$value.TSGs <-
          c(nTSGs:1, rep(0, length(tuson_TSGs$Gene) - nTSGs))
        
        colnames(tuson_OGs)[1] <- "gene"
        colnames(tuson_TSGs)[1] <- "gene"
        
        genes_in_bin <-
          left_join(genes_in_bin, tuson_OGs[, c(1, 44)], by = "gene")
        genes_in_bin <-
          left_join(genes_in_bin, tuson_TSGs[, c(1, 44)], by = "gene")
      }
      
      
      
      if (tumor_type == "BRCA") {
        genes_in_bin <- genes_in_bin[!is.na(genes_in_bin$chr), ]
        
        # https://www.sciencedirect.com/science/article/pii/S0092867418302149?via%3Dihub#app2
        # table S2: Analysis of proliferation screens in HMEC (breast cancer cell line) and HPNE (pancreatic adenocarcinoma cell line)
        library1 <- read_tsv(file = "data/misc/Sack2018_Library1.tsv")
        library2 <- read_tsv(file = "data/misc/Sack2018_Library2.tsv")
        
        
        table(library1$HMEC.p.val <= 0.1 & library1$HMEC.Log2 > 1)
        table(library2$HMEC.p.val <= 0.1 & library2$HMEC.Log2 > 1)
        
        library1 <- library1[!is.na(library1$HMEC.Log2),]
        library2 <- library2[!is.na(library2$HMEC.Log2),]
        
        
        library1$GO.status <-
          ifelse(library1$HMEC.p.val <= 0.03 &
                   library1$HMEC.Log2 > 1,
                 "GO",
                 "noGO")
        colnames(library1)[2] <- "gene"
        library2$GO.status <-
          ifelse(library2$HMEC.p.val <= 0.03 &
                   library2$HMEC.Log2 > 1,
                 "GO",
                 "noGO")
        colnames(library2)[2] <- "gene"
        
        nGO <-
          length(unique(c(library1[library1$GO.status == "GO",]$gene, library2[library2$GO.status == "GO",]$gene)))
        GO <-
          unique(c(library1[library1$GO.status == "GO",]$gene, library2[library2$GO.status == "GO",]$gene))
        genes_in_bin$GO.status <-
          ifelse(genes_in_bin$gene %in% GO, "GO", "noGO")
        
        
        GO.scores <- data.frame()
        for (bin_id_merged36 in levels(factor(genes_in_bin$bin_id_merged36))) {
          go.score <-
            dim(genes_in_bin[genes_in_bin$bin_id_merged36 == bin_id_merged36 &
                               genes_in_bin$GO.status == "GO",])[1] / dim(genes_in_bin[genes_in_bin$bin_id_merged36 == bin_id_merged36,])[1]
          
          GO.scores <-
            rbind(GO.scores,
                  cbind(bin_id_merged36 = bin_id_merged36,
                        go.score = go.score))
        }
      } else if (tumor_type == "PAAD") {
        genes_in_bin <- genes_in_bin[!is.na(genes_in_bin$chr),]
        library1 <- read_tsv(file = "data/misc/Sack2018_Library1.tsv")
        library2 <- read_tsv(file = "data/misc/Sack2018_Library2.tsv")
        
        # table(library1$HPNE.FDR <= 0.01 & library1$HPNE.Log2 > 1)
        # table(library2$HPNE.FDR <= 0.01 & library2$HPNE.Log2 > 1)
        
        library1 <- library1[!is.na(library1$HPNE.Log2),]
        library2 <- library2[!is.na(library2$HPNE.Log2),]
        
        
        library1$GO.status <-
          ifelse(library1$HPNE.FDR <= 0.01 &
                   library1$HPNE.Log2 > 1,
                 "GO",
                 "noGO")
        colnames(library1)[2] <- "gene"
        library2$GO.status <-
          ifelse(library2$HPNE.FDR <= 0.01 &
                   library2$HPNE.Log2 > 1,
                 "GO",
                 "noGO")
        colnames(library2)[2] <- "gene"
        
        nGO <-
          length(unique(c(library1[library1$GO.status == "GO",]$gene, library2[library2$GO.status == "GO",]$gene)))
        GO <-
          unique(c(library1[library1$GO.status == "GO",]$gene, library2[library2$GO.status == "GO",]$gene))
        genes_in_bin$GO.status <-
          ifelse(genes_in_bin$gene %in% GO, "GO", "noGO")
        
        
        GO.scores <- data.frame()
        for (bin_id_merged36 in levels(factor(genes_in_bin$bin_id_merged36))) {
          go.score <-
            dim(genes_in_bin[genes_in_bin$bin_id_merged36 == bin_id_merged36 &
                               genes_in_bin$GO.status == "GO",])[1] / dim(genes_in_bin[genes_in_bin$bin_id_merged36 == bin_id_merged36,])[1]
          
          GO.scores <-
            rbind(GO.scores,
                  cbind(bin_id_merged36 = bin_id_merged36,
                        go.score = go.score))
        }
      }
      
      merged <-
        full_join(bin_values, GO.scores, by = "bin_id_merged36")
      
      
      OG.scores <- data.frame()
      for (bin_id_merged36 in levels(factor(genes_in_bin$bin_id_merged36))) {
        og.score <-
          dim(genes_in_bin[genes_in_bin$bin_id_merged36 == bin_id_merged36 &
                             genes_in_bin$OG.status == "OG",])[1] / dim(genes_in_bin[genes_in_bin$bin_id_merged36 == bin_id_merged36,])[1]
        og.score.pancancer <-
          dim(genes_in_bin[genes_in_bin$bin_id_merged36 == bin_id_merged36 &
                             genes_in_bin$OG.status.pancancer == "OG",])[1] / dim(genes_in_bin[genes_in_bin$bin_id_merged36 == bin_id_merged36,])[1]
        tsg.score <-
          dim(genes_in_bin[genes_in_bin$bin_id_merged36 == bin_id_merged36 &
                             genes_in_bin$TSG.status == "TSG",])[1] / dim(genes_in_bin[genes_in_bin$bin_id_merged36 == bin_id_merged36,])[1]
        tsg.score.pancancer <-
          dim(genes_in_bin[genes_in_bin$bin_id_merged36 == bin_id_merged36 &
                             genes_in_bin$TSG.status.pancancer == "TSG",])[1] / dim(genes_in_bin[genes_in_bin$bin_id_merged36 == bin_id_merged36,])[1]
        
        OG.scores <-
          rbind(
            OG.scores,
            cbind(
              bin_id_merged36 = bin_id_merged36,
              og.score = og.score,
              og.score.pancancer = og.score.pancancer,
              tsg.score = tsg.score,
              tsg.score.pancancer = tsg.score.pancancer
            )
          )
      }
      
      merged <- full_join(merged, OG.scores, by = "bin_id_merged36")
      
      merged$cna_freq_ampl <-
        ifelse(is.na(merged$cna_freq_ampl), 0, merged$cna_freq_ampl)
      
      if (CTR == T) {
        # CONTROL: swap amplifications
        if (tumor_type == "BRCA") {
          tumor_type_control <- "PAAD"
        } else if (tumor_type == "PAAD" | tumor_type == "COADREAD") {
          tumor_type_control <- "BRCA"
        }
        cat(
          "\n > CONTROL: Amplification frequency of",
          tumor_type,
          "is replace with",
          tumor_type_control,
          "\n\n"
        )
        bin_values <-
          read.table(
            file = paste0(
              "results/tables/02_produceStatistics/amplifications_",
              tumor_type_control,
              "_",
              merged_1,
              "Mbp_table.txt"
            )
          )
        bin_values$bin_id_merged36 <-
          paste0(bin_values$chr, "_", bin_values$resize)
        
        merged$cna_freq_ampl_CTR <- bin_values$cna_freq_ampl
        
        merged$og.go.score <-
          as.numeric(merged$og.score) + as.numeric(merged$go.score)
        merged$mutation.score <-
          1 - log10(merged$mutations_norm) / min(log10(merged$mutations_norm))
        merged$og.go.score.mutations.norm <-
          merged$mutation.score + (as.numeric(merged$og.go.score))
        
        test <- cor.test(merged$og.go.score.mutations.norm,
                         merged$cna_freq_ampl,
                         method = "spearman")
        ctr <- cor.test(merged$og.go.score.mutations.norm,
                        merged$cna_freq_ampl_CTR,
                        method = "spearman")
        # plot(merged$og.go.score.mutations.norm, merged$cna_freq_ampl)
        
        stat_estimates_CTR <- rbind(stat_estimates_CTR,
                                    cbind(
                                      tumor_type,
                                      rbind(test$estimate,
                                            ctr$estimate),
                                      condition = rbind("test",
                                                        "CTR_swap.amplification")
                                    ))
        stat_signif_CTR <- rbind(stat_signif_CTR,
                                 cbind(
                                   tumor_type,
                                   rbind(test$p.value,
                                         ctr$p.value),
                                   condition = rbind("test.p",
                                                     "CTR_swap.amplification.p")
                                 ))
      }
      
      # SPEARMAN ----
      cor.test(merged$cna_freq_ampl,
               as.numeric(merged$og.score),
               method = "spearman")
      cor.test(merged$cna_freq_ampl,
               as.numeric(merged$go.score),
               method = "spearman")
      mu_score <-
        cor.test(merged$cna_freq_ampl, log10(as.numeric(merged$mutations_norm)), method = "spearman")
      
      merged$og.go.score <-
        as.numeric(merged$og.score) + as.numeric(merged$go.score)
      og.go_score <- cor.test(merged$cna_freq_ampl,
                              as.numeric(merged$og.go.score),
                              method = "spearman")
      # plot(merged$cna_freq_ampl, as.numeric(merged$og.go.score))
      
      merged$mutation.score <-
        1 - log10(merged$mutations_norm) / min(log10(merged$mutations_norm))
      cor.test(merged$cna_freq_ampl, merged$mutation.score, method = "spearman")
      # plot(merged$cna_freq_ampl, as.numeric(merged$mutation.score))
      
      
      merged$og.go.score.mutations.norm <-
        merged$mutation.score + (as.numeric(merged$og.go.score))
      mu.og.go_score <- cor.test(merged$og.go.score.mutations.norm,
                                 merged$cna_freq_ampl,
                                 method = "spearman")
      # plot(merged$og.go.score.mutations.norm, merged$cna_freq_ampl)
      
      merged[, 14:22] <- apply(merged[, 14:22], 2, as.numeric)
      
      a1 <- ggscatter(
        merged,
        x = "og.score",
        y = "cna_freq_ampl",
        add = "reg.line",
        conf.int = T,
        cor.coef = T,
        cor.method = "spearman",
        # size = "og.score", color = "mutation.score",
        xlab = "(OG)score",
        ylab = "Amplification frequency",
        cor.coeff.args = list(method = "spearman")
      )  +
        ggtitle(label = tumor_type,
                subtitle = paste("ampl.f ~ OG.score"))
      
      a2 <- ggscatter(
        merged,
        x = "go.score",
        y = "cna_freq_ampl",
        add = "reg.line",
        conf.int = T,
        cor.coef = T,
        cor.method = "spearman",
        # size = "og.score", color = "mutation.score",
        xlab = "(GO)score",
        ylab = "Amplification frequency",
        cor.coeff.args = list(method = "spearman")
      )  +
        ggtitle(label = tumor_type,
                subtitle = paste("ampl.f ~ GO.score"))
      
      a3 <- ggscatter(
        merged,
        x = "og.go.score",
        y = "cna_freq_ampl",
        add = "reg.line",
        conf.int = T,
        cor.coef = T,
        cor.method = "spearman",
        # size = "og.score", color = "mutation.score",
        xlab = "(OG+GO)score",
        ylab = "Amplification frequency",
        cor.coeff.args = list(method = "spearman")
      )   +
        ggtitle(label = tumor_type,
                subtitle = paste("ampl.f ~ OG+GO.score"))
      
      a4 <- ggscatter(
        merged,
        x = "mutation.score",
        y = "cna_freq_ampl",
        add = "reg.line",
        conf.int = T,
        cor.coef = T,
        cor.method = "spearman",
        # size = "og.score", color = "mutation.score",
        xlab = "MUT.score",
        ylab = "Amplification frequency",
        cor.coeff.args = list(method = "spearman")
      )   +
        ggtitle(label = tumor_type,
                subtitle = paste("ampl.f ~ MUT.score"))
      
      (
        a5 <-
          ggscatter(
            merged,
            x = "og.go.score.mutations.norm",
            y = "cna_freq_ampl",
            add = "reg.line",
            conf.int = T,
            cor.coef = T,
            cor.method = "spearman",
            # size = "mutation.score",
            # color = "og.score",
            xlab = "(OG+GO).score*MUT.score",
            ylab = "Amplification frequency",
            cor.coeff.args = list(method = "spearman")
          ) +
          ggtitle(
            label = paste0(tumor_type, " - most complete model"),
            subtitle = paste("ampl.f ~ (OG+GO).score*MUT.score")
          )
      )
      
      a6 <- ggscatter(
        merged,
        x = "og.go.score",
        y = "mutation.score",
        add = "reg.line",
        conf.int = T,
        cor.coef = T,
        cor.method = "spearman",
        xlab = "(OG+GO).score",
        ylab = "MUT.score",
        cor.coeff.args = list(method = "spearman")
      ) +
        ggtitle(label = tumor_type,
                subtitle = paste("CONTROL: MUT.score ~ (OG+GO).score"))
      a7 <- ggscatter(
        merged,
        x = "go.score",
        y = "mutation.score",
        add = "reg.line",
        conf.int = T,
        cor.coef = T,
        cor.method = "spearman",
        xlab = "(GO).score",
        ylab = "MUT.score",
        cor.coeff.args = list(method = "spearman")
      ) +
        ggtitle(label = tumor_type,
                subtitle = paste("CONTROL: MUT.score ~ (GO).score"))
      a8 <- ggscatter(
        merged,
        x = "og.score",
        y = "mutation.score",
        add = "reg.line",
        conf.int = T,
        cor.coef = T,
        cor.method = "spearman",
        xlab = "(OG).score",
        ylab = "MUT.score",
        cor.coeff.args = list(method = "spearman")
      ) +
        ggtitle(label = tumor_type,
                subtitle = paste("CONTROL: MUT.score ~ (OG).score"))
      
      pdf(
        file = paste0(
          "results/plots/000_paper_plots/00_SupplementaryFig3_",
          tumor_type,
          ".pdf"
        ),
        width = 12,
        height = 7
      )
      print(ggarrange(a1, a2))
      print(ggarrange(a3, a4))
      print(a5)
      print(ggarrange(a6, a7, a8, ncol = 3))
      dev.off()
      
      stat_estimates <- rbind(stat_estimates,
                              cbind(
                                tumor_type,
                                rbind(
                                  og.go_score$estimate,
                                  mu_score$estimate,
                                  mu.og.go_score$estimate
                                ),
                                condition = rbind("og.go_corS",
                                                  "mu_corS",
                                                  "mu.og.go_corS")
                              ))
      stat_signif <- rbind(stat_signif,
                           cbind(
                             tumor_type,
                             rbind(
                               og.go_score$p.value,
                               mu_score$p.value,
                               mu.og.go_score$p.value
                             ),
                             condition = rbind("og.go_p.corS",
                                               "mu_p.corS",
                                               "mu.og.go_pcorS")
                           ))
    }
  })
})

colnames(stat_estimates) <- c("tumor_types", "rho", "condition")
stat_estimates$rho <- as.numeric(stat_estimates$rho)

colnames(stat_signif) <- c("tumor_types", "p.val", "condition")
stat_signif$p.val <- as.numeric(stat_signif$p.val)

p1 <-
  ggplot(stat_estimates, aes(x = as.factor(tumor_types), y = rho)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(condition))) +
  theme_classic() + scale_fill_brewer(palette = "Blues")

p2 <-
  ggplot(stat_signif, aes(x = as.factor(tumor_types), y = log10(p.val))) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(condition))) +
  theme_classic() + scale_fill_brewer(palette = "Reds")

cat(" \n > Producing Fig. 3 \n\n")
pdf("results/plots/000_paper_plots/00_Fig3.pdf", height = 10)
print(ggarrange(p1, p2, ncol = 1))
dev.off()

if (CTR == T) {
  cat(" \n > Producing Fig. S4 \n\n")
  
  colnames(stat_estimates_CTR) <-
    c("tumor_types", "rho", "condition")
  stat_estimates_CTR$rho <- as.numeric(stat_estimates_CTR$rho)
  
  colnames(stat_signif_CTR) <-
    c("tumor_types", "p.val", "condition")
  stat_signif_CTR$p.val <- as.numeric(stat_signif_CTR$p.val)
  
  p1 <-
    ggplot(stat_estimates_CTR, aes(x = as.factor(tumor_types), y = rho)) +
    geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(condition))) +
    theme_classic() + scale_fill_brewer(palette = "Blues")
  
  p2 <-
    ggplot(stat_signif_CTR, aes(x = as.factor(tumor_types), y = log10(p.val))) +
    geom_bar(stat = "identity", position = "dodge", aes(fill = as.factor(condition))) +
    theme_classic() + scale_fill_brewer(palette = "Reds")
  
  pdf("results/plots/000_paper_plots/00_SupplementaryFig4.pdf",
      height = 10)
  print(ggarrange(p1, p2, ncol = 1))
  dev.off()
}

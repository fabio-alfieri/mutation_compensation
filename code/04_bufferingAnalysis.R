# BUFFERING ANALYSIS

suppressMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(dplyr)
  library(parallel)
})

cat("\n\n > This script provides the degree of buffering of mutations within different segment's allele copy number \n\n")

setwd("../")
# setwd("/home/ieo5099/Desktop_linux/mutation_compensation/")

mutations_wCCF <-
  fread(paste0("data/misc/PCAWG_mutations_wCCF_allTumorTypes.tsv"), sep = "\t",
        nThread = 8)

mutations_wCCF$multiplicity <-
  as.factor(mutations_wCCF$multiplicity)

# define coding and non-coding mutations ----
mutations_wCCF$multiplicity <-
  as.numeric(as.character(mutations_wCCF$multiplicity))

coding <- mutations_wCCF[mutations_wCCF$consequence_type %in%
                           names(table(mutations_wCCF$consequence_type)[str_detect(names(table(mutations_wCCF$consequence_type)),
                                                                                   "exon")]) |
                           mutations_wCCF$consequence_type %in%
                           names(table(mutations_wCCF$consequence_type)[str_detect(names(table(mutations_wCCF$consequence_type)),
                                                                                   "missense")]) |
                           mutations_wCCF$consequence_type %in%
                           names(table(mutations_wCCF$consequence_type)[str_detect(names(table(mutations_wCCF$consequence_type)),
                                                                                   "synonymous")]),]

noncoding <- mutations_wCCF[mutations_wCCF$consequence_type  %in%
                              names(table(mutations_wCCF$consequence_type)[str_detect(names(table(mutations_wCCF$consequence_type)),
                                                                                      "intergenic")]) |
                              mutations_wCCF$consequence_type %in%
                              names(table(mutations_wCCF$consequence_type)[str_detect(names(table(mutations_wCCF$consequence_type)),
                                                                                      "intron")]) & 
                              !(mutations_wCCF$consequence_type %in%
                                  names(table(mutations_wCCF$consequence_type)[str_detect(names(table(mutations_wCCF$consequence_type)),
                                                                                          "exon")]) |
                                  mutations_wCCF$consequence_type %in%
                                  names(table(mutations_wCCF$consequence_type)[str_detect(names(table(mutations_wCCF$consequence_type)),
                                                                                          "missense")]) |
                                  mutations_wCCF$consequence_type %in%
                                  names(table(mutations_wCCF$consequence_type)[str_detect(names(table(mutations_wCCF$consequence_type)),
                                                                                          "synonymous")])),]

timing_kar <- function(ta) {
  ta1.0 <- ta[ta$karyotype == "1:0", ]
  ta1.0$timing <- ifelse(ta1.0$multiplicity <= 2, "ns", "ns")
  ta1.1 <- ta[ta$karyotype == "1:1", ]
  ta1.1$timing <- ifelse(ta1.1$multiplicity <= 2, "ns", "ns")
  
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
  
  ta <-
    rbind(
      ta1.0,
      ta1.1,
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

buffering_wdiploid_kar <- function(ta) {
  ta1.0 <- ta[ta$karyotype == "1:0", ]
  ta1.0$buffering <- ifelse(ta1.0$multiplicity <= 2, "deleted_unbuffered", "deleted_unbuffered")
  ta1.1 <- ta[ta$karyotype == "1:1", ]
  ta1.1$buffering <- ifelse(ta1.1$multiplicity <= 2, "diploid_unbuffered", "diploid_unbuffered")
  
  ta2.0 <- ta[ta$karyotype == "2:0", ]
  ta2.0$buffering <- ifelse(ta2.0$multiplicity <= 2, "cn_loh_unbuffered", "cn_loh_unbuffered")
  ta2.1 <- ta[ta$karyotype == "2:1", ]
  ta2.1$buffering <- ifelse(ta2.1$multiplicity == 1, "amplified_buffered", "amplified_unbuffered")
  ta2.2 <- ta[ta$karyotype == "2:2", ]
  ta2.2$buffering <- ifelse(ta2.2$multiplicity == 1, "amplified_buffered", "amplified_buffered")
  
  ta3.0 <- ta[ta$karyotype == "3:0", ]
  ta3.0$buffering <- ifelse(ta3.0$multiplicity == 1, "amplified_buffered", "amplified_unbuffered")
  ta3.1 <- ta[ta$karyotype == "3:1", ]
  ta3.1$buffering <- ifelse(ta3.1$multiplicity <= 2, "amplified_buffered", "amplified_unbuffered")
  ta3.2 <- ta[ta$karyotype == "3:2", ]
  ta3.2$buffering <- ifelse(ta3.2$multiplicity <= 2, "amplified_buffered", "amplified_buffered")
  ta3.3 <- ta[ta$karyotype == "3:3", ]
  ta3.3$buffering <- ifelse(ta3.3$multiplicity <= 2, "amplified_buffered", "amplified_buffered")
  
  ta4.0 <- ta[ta$karyotype == "4:0", ]
  ta4.0$buffering <- ifelse(ta4.0$multiplicity <= 2, "amplified_buffered", "amplified_unbuffered")
  ta4.1 <- ta[ta$karyotype == "4:1", ]
  ta4.1$buffering <- ifelse(ta4.1$multiplicity <= 3, "amplified_buffered", "amplified_unbuffered")
  ta4.2 <- ta[ta$karyotype == "4:2", ]
  ta4.2$buffering <- ifelse(ta4.2$multiplicity <= 3, "amplified_buffered", "amplified_buffered")
  ta4.3 <- ta[ta$karyotype == "4:3", ]
  ta4.3$buffering <- ifelse(ta4.3$multiplicity <= 3, "amplified_buffered", "amplified_buffered")
  
  ta5.0 <- ta[ta$karyotype == "5:0", ]
  ta5.0$buffering <- ifelse(ta5.0$multiplicity <= 3, "amplified_buffered", "amplified_unbuffered")
  ta5.1 <- ta[ta$karyotype == "5:1", ]
  ta5.1$buffering <- ifelse(ta5.1$multiplicity <= 4, "amplified_buffered", "amplified_unbuffered")
  ta5.2 <- ta[ta$karyotype == "5:2", ]
  ta5.2$buffering <- ifelse(ta5.2$multiplicity <= 4, "amplified_buffered", "amplified_buffered")
  
  ta <-
    rbind(
      ta1.0,
      ta1.1,
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

buffering_kar <- function(ta) {
  ta1.0 <- ta[ta$karyotype == "1:0", ]
  ta1.0$buffering <- ifelse(ta1.0$multiplicity <= 2, "unbuffered", "unbuffered")
  ta1.1 <- ta[ta$karyotype == "1:1", ]
  ta1.1$buffering <- ifelse(ta1.1$multiplicity <= 2, "unbuffered", "unbuffered")
  
  ta2.0 <- ta[ta$karyotype == "2:0", ]
  ta2.0$buffering <- ifelse(ta2.0$multiplicity <= 2, "unbuffered", "unbuffered")
  ta2.1 <- ta[ta$karyotype == "2:1", ]
  ta2.1$buffering <- ifelse(ta2.1$multiplicity == 1, "buffered", "unbuffered")
  ta2.2 <- ta[ta$karyotype == "2:2", ]
  ta2.2$buffering <- ifelse(ta2.2$multiplicity == 1, "buffered", "buffered")
  
  ta3.0 <- ta[ta$karyotype == "3:0", ]
  ta3.0$buffering <- ifelse(ta3.0$multiplicity == 1, "buffered", "unbuffered")
  ta3.1 <- ta[ta$karyotype == "3:1", ]
  ta3.1$buffering <- ifelse(ta3.1$multiplicity <= 2, "buffered", "unbuffered")
  ta3.2 <- ta[ta$karyotype == "3:2", ]
  ta3.2$buffering <- ifelse(ta3.2$multiplicity <= 2, "buffered", "buffered")
  ta3.3 <- ta[ta$karyotype == "3:3", ]
  ta3.3$buffering <- ifelse(ta3.3$multiplicity <= 2, "buffered", "buffered")
  
  ta4.0 <- ta[ta$karyotype == "4:0", ]
  ta4.0$buffering <- ifelse(ta4.0$multiplicity <= 2, "buffered", "unbuffered")
  ta4.1 <- ta[ta$karyotype == "4:1", ]
  ta4.1$buffering <- ifelse(ta4.1$multiplicity <= 3, "buffered", "unbuffered")
  ta4.2 <- ta[ta$karyotype == "4:2", ]
  ta4.2$buffering <- ifelse(ta4.2$multiplicity <= 3, "buffered", "buffered")
  ta4.3 <- ta[ta$karyotype == "4:3", ]
  ta4.3$buffering <- ifelse(ta4.3$multiplicity <= 3, "buffered", "buffered")
  
  ta5.0 <- ta[ta$karyotype == "5:0", ]
  ta5.0$buffering <- ifelse(ta5.0$multiplicity <= 3, "buffered", "unbuffered")
  ta5.1 <- ta[ta$karyotype == "5:1", ]
  ta5.1$buffering <- ifelse(ta5.1$multiplicity <= 4, "buffered", "unbuffered")
  ta5.2 <- ta[ta$karyotype == "5:2", ]
  ta5.2$buffering <- ifelse(ta5.2$multiplicity <= 4, "buffered", "buffered")
  
  ta <-
    rbind(
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

coding_timing <- buffering_kar(coding)

coding_timing <- coding_timing[!is.na(coding_timing$buffering),]

final_coding <- full_join(as.data.frame(table(coding_timing$project_code)),
                          as.data.frame(table(coding_timing$project_code, coding_timing$buffering)), by = "Var1")

final_coding$prop <- final_coding$Freq.y/final_coding$Freq.x

ggplot(final_coding, aes(x = Var1, y = prop, fill = Var2)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "top") +
  ggtitle("Coding mutations") +
  theme_classic()


noncoding_timing <- buffering_kar(noncoding)
noncoding_timing <- noncoding_timing[!is.na(noncoding_timing$buffering),]

ft <- fisher.test(rbind(table(coding_timing$buffering), table(noncoding_timing$buffering)))
ft$p.value

final_noncoding <- full_join(as.data.frame(table(noncoding_timing$project_code)),
                             as.data.frame(table(noncoding_timing$project_code, 
                                                 noncoding_timing$buffering)), by = "Var1")

final_noncoding$prop <- final_noncoding$Freq.y/final_noncoding$Freq.x

ggplot(final_noncoding, aes(x = Var1, y = prop, fill = Var2)) +
  geom_bar(stat = "identity") + 
  theme(legend.position = "top") +
  ggtitle("Non-coding mutations")


final_merged <- rbind(cbind(final_noncoding, status = "noncoding"), 
                      cbind(final_coding, status = "coding"))
final_merged$xlabel <- paste0(final_merged$Var1,"_",final_merged$status)

library(qualpalr)
set.seed(2)
pal = qualpal(length(levels(factor(final_merged$Var2)))+12, colorspace = "pretty")

ggplot(final_merged, aes(x = Var1, y = prop/2, fill = Var2)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_fill_manual(values=sample(pal$hex, 2)) +
  ggtitle("PCAWG: coding and non-coding mutations in all regions (amplified, deleted and diploid)") +
  xlab("Tumor types") +
  ylab("%")


# still wo diploid ----
noncoding_timing <- timing_kar(noncoding_timing)
coding_timing <- timing_kar(coding_timing)

coding_timing <- coding_timing[!is.na(coding_timing$timing),]
noncoding_timing <- noncoding_timing[!is.na(noncoding_timing$timing),]

final_coding <- full_join(as.data.frame(table(coding_timing$project_code)),
                          as.data.frame(table(coding_timing$project_code, paste0(coding_timing$timing, "_", coding_timing$buffering))), by = "Var1")
final_noncoding <- full_join(as.data.frame(table(noncoding_timing$project_code)),
                             as.data.frame(table(noncoding_timing$project_code, 
                                                 paste0(noncoding_timing$timing,"_",noncoding_timing$buffering))), by = "Var1")

final_merged_timing <- rbind(cbind(final_noncoding, status = "noncoding"), 
                             cbind(final_coding, status = "coding"))
final_merged_timing$xlabel <- paste0(final_merged_timing$Var1,"_",final_merged_timing$status)

final_merged_timing$prop <- final_merged_timing$Freq.y/final_merged_timing$Freq.x

library(qualpalr)
set.seed(4)
pal = qualpal(length(levels(factor(final_merged_timing$Var2)))+10, colorspace = "pretty")
# colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

pdf(file = "results/plots/04_bufferingAnalysis/buffered_unbuffered_later_early.pdf", width = 10, height = 8)
ggplot(final_merged_timing, aes(x = Var1, y = prop/2, fill = Var2)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_fill_manual(values=sample(pal$hex), 4) +
  ggtitle("PCAWG: coding and non-coding mutations in amplified regions") +
  xlab("Tumor types") +
  ylab("%")
dev.off()

# w diploid ----
noncoding_timing <- buffering_wdiploid_kar(noncoding)
coding_timing <- buffering_wdiploid_kar(coding)

noncoding_timing <- timing_kar(noncoding_timing)
coding_timing <- timing_kar(coding_timing)

coding_timing <- coding_timing[!is.na(coding_timing$timing),]
noncoding_timing <- noncoding_timing[!is.na(noncoding_timing$timing),]
coding_timing <- coding_timing[!is.na(coding_timing$buffering),]
noncoding_timing <- noncoding_timing[!is.na(noncoding_timing$buffering),]

merged <- rbind(cbind(noncoding_timing, status = "noncoding"), 
                cbind(coding_timing, status = "coding"))


final_coding <- full_join(as.data.frame(table(coding_timing$project_code)),
                          as.data.frame(table(coding_timing$project_code, paste0(coding_timing$timing, "_", coding_timing$buffering))), by = "Var1")
final_noncoding <- full_join(as.data.frame(table(noncoding_timing$project_code)),
                             as.data.frame(table(noncoding_timing$project_code, 
                                                 paste0(noncoding_timing$timing,"_",noncoding_timing$buffering))), by = "Var1")

final_merged_timing <- rbind(cbind(final_noncoding, status = "noncoding"), 
                             cbind(final_coding, status = "coding"))
final_merged_timing$xlabel <- paste0(final_merged_timing$Var1,"_",final_merged_timing$status)

final_merged_timing$prop <- final_merged_timing$Freq.y/final_merged_timing$Freq.x

library(qualpalr)
set.seed(1)
pal = qualpal(length(levels(factor(final_merged_timing$Var2))), colorspace = "pretty")
# colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

pdf(file = "results/plots/04_bufferingAnalysis/buffered_unbuffered_later_early_w_diploid.pdf", width = 10, height = 8)
ggplot(final_merged_timing, aes(x = Var1, y = prop/2, fill = Var2)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_fill_manual(values=pal$hex) +
  ggtitle("PCAWG: coding and non-coding mutations in all regions (amplified, deleted and diploid)") +
  xlab("Tumor types") +
  ylab("%")
dev.off()
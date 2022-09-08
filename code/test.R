#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(c("-a", "--analysis"), type="character", default=NULL, 
              help="Options are:
              \t - 00 or geneLevelAnalysis
              \t - 01 or binLevelAnalysis
              \t - 02 or produceStatistics
              \t - 03 or BRCAandPAADAnalysis
              \t - 04 or timingAnalysis
              \t - 05 or protectionAnalysis
              \t - all", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$analysis)) {
  print_help(opt_parser)
  stop("please specify the analysis you want to perform!", call.=FALSE)
}
if (!any(opt$analysis %in% c("00","geneLevelAnalysis",
                             "01","binLevelAnalysis",
                             "02","statisticalAnalysis",
                             "03","BRCAandPAADAnalysis",
                             "04","timingAnalysis",
                             "05","protectionAnalysis",
                             "all"
                             ))) {
  print_help(opt_parser)
  stop("typo in the analysis flag, plase see above for the available options!", call.=FALSE)
}

if (opt$analysis == "geneLevelAnalysis" | opt$analysis == "00") {
  system(paste0("Rscript --vanilla code/00_geneLevelAnalysis.R"))
}

if (opt$analysis == "binLevelAnalysis" | opt$analysis == "01") {
  system(paste0("Rscript --vanilla code/01_binLevelAnalysis.R"))
}

if (opt$analysis == "protectionAnalysis" | opt$analysis == "05") {
  system(paste0("Rscript --vanilla code/05_protectionAnalysis.R"))
  stop()
}

# if (opt$analysis == "BRCAandPAADAnalysis" | opt$analysis == "03") {
#   system(paste0("Rscript --vanilla code/02_statisticalAnalysis.R"))
# }
# 
# if (opt$analysis == "statisticalAnalysis" | opt$analysis == "02") {
#   system(paste0("Rscript --vanilla code/02_statisticalAnalysis.R"))
# }
# 
# if (opt$analysis == "statisticalAnalysis" | opt$analysis == "02") {
#   system(paste0("Rscript --vanilla code/02_statisticalAnalysis.R"))
# }5h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[1mRows: [22m[34m1912[39m [1mColumns: [22m[34m1[39m
[36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────────[39m
[1mDelimiter:[22m ","
[31mchr[39m (1): gene

[36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
[36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.
[?25h[1mRows: [22m[34m781[39m [1mColumns: [22m[34m1[39m
[36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────────────────────────────────────────────────[39m
[1mDelimiter:[22m ","
[31mchr[39m (1): gene

[36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
[36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.
[?25h
FALSE  TRUE 
 2115   328 
[?25h
FALSE  TRUE 
 2401    42 
[?25h

 >> Contingency tables 

[?25h     [,1] [,2]
TRUE  328   42
TRUE   46  120
[?25h[?25h     FALSE TRUE
[1,]  2401   42
[2,]   914  120
[?25h
	Fisher's Exact Test for Count Data

data:  rbind(table(protected$V1 %in% nonessential$gene), table(unprotected$V1 %in% nonessential$gene))
p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  5.189696 11.021757
sample estimates:
odds ratio 
  7.499359 

[?25h     FALSE TRUE
[1,]  2115  328
[2,]   988   46
[?25h
	Fisher's Exact Test for Count Data

data:  rbind(table(protected$V1 %in% common_essential$gene), table(unprotected$V1 %in% common_essential$gene))
p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.2136024 0.4140493
sample estimates:
odds ratio 
 0.3003023 

[?25h[?25hError: 
Execution halted
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ^C^C
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ^C
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ 
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ^C
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ^C
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ 
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ Rscript test.R -a 05[K[K[K[K[K[K[K[K[K[K[K[K
code/        data/        data.zip     results/     results.zip  Rplots.pdf   test.R       
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ Rscript code/
00_geneLevelAnalysis.R    02_statisticalAnalysis.R  04_timingAnalysis.R       .Rapp.history
01_binLevelAnalysis.R     03_BRCA-PAAD_analysis.R   05_protectionAnalysis.R   
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ Rscript code/05_protectionAnalysis.R 


 >> Loading libraries 

^C
Execution halted
[?25h(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ..
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  2 12:04  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 11:59  [34;42mmutation_compensation[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll..Rscript code/05_protectionAnalysis.R 
Fatal error: cannot open file 'code/05_protectionAnalysis.R': No such file or directory
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ cd mutation_compensation/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ cd[K[Kll
total 6,6G
drwxrwxrwx 1 laninst laninst 4,0K set  7 11:59 [0m[34;42m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  2 12:04 [01;34m..[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  6 17:59 [01;34mcode[0m/
drwxrwxr-x 1 laninst laninst 4,0K ago 31 11:52 [01;34mdata[0m/
-rw-rw-r-- 1 laninst laninst 5,4G set  1 16:38 [01;31mdata.zip[0m
drwxrwxr-x 1 laninst laninst 4,0K ago 23 10:48 [01;34mresults[0m/
-rw-rw-r-- 1 laninst laninst 1,3G set  1 16:25 [01;31mresults.zip[0m
-rw-rw-r-- 1 laninst laninst  24K set  7 17:32 Rplots.pdf
-rw-r--r-- 1 laninst laninst 6,1K set  7 17:34 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ cd code/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ ll
total 128K
drwxr-xr-x 1 laninst laninst 4,0K set  6 17:59 [0m[01;34m.[0m/
drwxrwxrwx 1 laninst laninst 4,0K set  7 11:59 [34;42m..[0m/
-rwxr-xr-x 1 laninst laninst 8,7K set  6 18:36 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:33 [01;32m05_protectionAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst    0 ago 30 15:03 [01;32m.Rapp.history[0m*
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ cd [K[K[KRscript --vai[Knilla 095[K[K0[K5_protectionAnalysis.R 


 >> Loading libraries 

[?25h[?25h

 >> Calculate protection index 

[?25h[?25h		 LUAD 
Error in file(file, "rt") : cannot open the connection
Calls: suppressMessages ... withCallingHandlers -> read.csv -> read.table -> file
In addition: Warning message:
In file(file, "rt") :
  cannot open file 'results/tables/00_geneLevelAnalysis/LUAD_geneLevel.csv': No such file or directory
Execution halted
[?25h(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ Rscript --vanilla 05_protectionAnalysis.R 


 >> Loading libraries 

[?25h[?25h[?25h

 >> Calculate protection index 

[?25h[?25h		 LUAD 
		 LUSC 
		 BRCA 
		 CESC 
		 THCA 
		 HNSC 
		 PAAD 
		 GBMLGG 
		 COADREAD 
[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h

 >> Extacting common protected and unprotected genes across tumor types 

[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h

 >> Calculating essentiality scores for each gene of the two categories 

[?25h[?25h[?25h[?25h
	Welch Two Sample t-test

data:  mean_crispr_effect[mean_crispr_effect$genes %in% protected$V1, ]$mean and mean_crispr_effect[mean_crispr_effect$genes %in% unprotected$V1, ]$mean
t = -9.6362, df = 2647.7, p-value < 2.2e-16
alternative hypothesis: true difference in means is less than 0
95 percent confidence interval:
        -Inf -0.09674922
sample estimates:
  mean of x   mean of y 
-0.18701814 -0.07034669 

[?25h[?25h[?25h[?25h
	Wilcoxon rank sum test with continuity correction

data:  mean_crispr_effect[mean_crispr_effect$genes %in% protected$V1, ]$mean and mean_crispr_effect[mean_crispr_effect$genes %in% unprotected$V1, ]$mean
W = 828944, p-value < 2.2e-16
alternative hypothesis: true location shift is less than 0

[?25h[?25h

 >> Gene Ontology analysis 

[?25h > Analyzing the PROTECTED gene set 

[?25h[?25h'select()' returned 1:many mapping between keys and columns
[?25h'select()' returned 1:1 mapping between keys and columns
[?25h'select()' returned 1:many mapping between keys and columns
[?25h[?25h           ONTOLOGY                                            Description
GO:0033036       BP                             macromolecule localization
GO:0008104       BP                                   protein localization
GO:0051641       BP                                  cellular localization
GO:0022613       BP                   ribonucleoprotein complex biogenesis
GO:0034613       BP                          cellular protein localization
GO:0070727       BP                    cellular macromolecule localization
GO:0007005       BP                             mitochondrion organization
GO:0051668       BP                           localization within membrane
GO:0045184       BP                  establishment of protein localization
GO:0043604       BP                             amide biosynthetic process
GO:0043603       BP                       cellular amide metabolic process
GO:0043933       BP                protein-containing complex organization
GO:0015031       BP                                      protein transport
GO:0042254       BP                                    ribosome biogenesis
GO:0050790       BP                       regulation of catalytic activity
GO:0046907       BP                                intracellular transport
GO:0051336       BP                       regulation of hydrolase activity
GO:0061024       BP                                  membrane organization
GO:0006518       BP                              peptide metabolic process
GO:0006364       BP                                        rRNA processing
GO:1901566       BP           organonitrogen compound biosynthetic process
GO:0002181       BP                                cytoplasmic translation
GO:0019637       BP                      organophosphate metabolic process
GO:0006412       BP                                            translation
GO:0009056       BP                                      catabolic process
GO:0051246       BP                regulation of protein metabolic process
GO:0043043       BP                           peptide biosynthetic process
GO:0072657       BP                       protein localization to membrane
GO:0065003       BP                    protein-containing complex assembly
GO:0071702       BP                            organic substance transport
GO:0010256       BP                       endomembrane system organization
GO:0016072       BP                                 rRNA metabolic process
GO:0044248       BP                             cellular catabolic process
GO:0055086       BP nucleobase-containing small molecule metabolic process
GO:0051649       BP                  establishment of localization in cell
GO:0034470       BP                                       ncRNA processing
GO:0052548       BP                   regulation of endopeptidase activity
GO:1901575       BP                    organic substance catabolic process
GO:1990778       BP                 protein localization to cell periphery
GO:0043085       BP              positive regulation of catalytic activity
GO:0044255       BP                       cellular lipid metabolic process
GO:0052547       BP                       regulation of peptidase activity
GO:0006644       BP                         phospholipid metabolic process
GO:0071705       BP                            nitrogen compound transport
GO:0006508       BP                                            proteolysis
GO:0072659       BP                protein localization to plasma membrane
GO:0055088       BP                                      lipid homeostasis
GO:0044281       BP                       small molecule metabolic process
GO:0051247       BP       positive regulation of protein metabolic process
GO:0006753       BP                 nucleoside phosphate metabolic process
GO:0033108       BP       mitochondrial respiratory chain complex assembly
GO:0006629       BP                                lipid metabolic process
GO:0030162       BP                              regulation of proteolysis
GO:0032268       BP       regulation of cellular protein metabolic process
GO:0009117       BP                           nucleotide metabolic process
GO:0006793       BP                           phosphorus metabolic process
GO:0008610       BP                             lipid biosynthetic process
GO:0006886       BP                        intracellular protein transport
GO:0006796       BP        phosphate-containing compound metabolic process
GO:0098876       BP      vesicle-mediated transport to the plasma membrane
GO:0006396       BP                                         RNA processing
GO:0044093       BP              positive regulation of molecular function
GO:0090407       BP                   organophosphate biosynthetic process
GO:0022904       BP                   respiratory electron transport chain
GO:0042157       BP                          lipoprotein metabolic process
GO:0005739       CC                                          mitochondrion
GO:0031967       CC                                     organelle envelope
GO:0031975       CC                                               envelope
GO:1902494       CC                                      catalytic complex
GO:0098798       CC               mitochondrial protein-containing complex
GO:0005740       CC                                 mitochondrial envelope
GO:0031966       CC                                 mitochondrial membrane
GO:0005840       CC                                               ribosome
GO:0005743       CC                           mitochondrial inner membrane
GO:1990904       CC                              ribonucleoprotein complex
GO:0005730       CC                                              nucleolus
GO:0098800       CC           inner mitochondrial membrane protein complex
GO:0019866       CC                               organelle inner membrane
GO:0044391       CC                                      ribosomal subunit
GO:0070062       CC                                  extracellular exosome
GO:0005635       CC                                       nuclear envelope
GO:0043230       CC                                extracellular organelle
GO:0065010       CC               extracellular membrane-bounded organelle
GO:1903561       CC                                  extracellular vesicle
GO:0140535       CC               intracellular protein-containing complex
GO:0005783       CC                                  endoplasmic reticulum
GO:1990234       CC                                    transferase complex
GO:0005759       CC                                   mitochondrial matrix
GO:0048471       CC                        perinuclear region of cytoplasm
GO:1990204       CC                                 oxidoreductase complex
GO:0098796       CC                               membrane protein complex
GO:0005746       CC                              mitochondrial respirasome
GO:0005815       CC                          microtubule organizing center
GO:0005856       CC                                           cytoskeleton
GO:0005777       CC                                             peroxisome
GO:0042579       CC                                              microbody
GO:0000313       CC                                    organellar ribosome
GO:0005761       CC                                 mitochondrial ribosome
GO:0098803       CC                              respiratory chain complex
GO:0015934       CC                                large ribosomal subunit
                 pvalue     p.adjust Count
GO:0033036 2.925955e-09 1.980579e-05   362
GO:0008104 9.178469e-08 2.533672e-04   315
GO:0051641 1.122916e-07 2.533672e-04   353
GO:0022613 1.656384e-07 2.803016e-04    70
GO:0034613 2.945955e-07 3.323528e-04   230
GO:0070727 2.945955e-07 3.323528e-04   230
GO:0007005 6.662274e-07 6.442419e-04   110
GO:0051668 1.389792e-06 1.175938e-03   107
GO:0045184 2.213014e-06 1.371617e-03   241
GO:0043604 2.260198e-06 1.371617e-03   135
GO:0043603 2.352840e-06 1.371617e-03   173
GO:0043933 2.523639e-06 1.371617e-03   234
GO:0015031 2.997987e-06 1.371617e-03   233
GO:0042254 3.184505e-06 1.371617e-03    55
GO:0050790 3.202132e-06 1.371617e-03   320
GO:0046907 3.242115e-06 1.371617e-03   220
GO:0051336 3.931395e-06 1.525623e-03   136
GO:0061024 4.056907e-06 1.525623e-03   114
GO:0006518 4.399554e-06 1.567399e-03   139
GO:0006364 9.264776e-06 3.135663e-03    40
GO:1901566 1.203711e-05 3.879961e-03   240
GO:0002181 1.768690e-05 5.237549e-03    44
GO:0019637 1.779637e-05 5.237549e-03   100
GO:0006412 2.037548e-05 5.513444e-03   114
GO:0009056 2.070448e-05 5.513444e-03   313
GO:0051246 2.117736e-05 5.513444e-03   330
GO:0043043 2.291918e-05 5.711135e-03   117
GO:0072657 2.362414e-05 5.711135e-03    95
GO:0065003 2.481370e-05 5.791861e-03   208
GO:0071702 4.545409e-05 1.025596e-02   313
GO:0010256 5.380128e-05 1.174777e-02    80
GO:0016072 6.540618e-05 1.342309e-02    45
GO:0044248 6.696129e-05 1.342309e-02   270
GO:0055086 6.742283e-05 1.342309e-02    63
GO:0051649 8.562687e-05 1.656024e-02   274
GO:0034470 8.983584e-05 1.689163e-02    66
GO:0052548 1.126131e-04 2.060210e-02    61
GO:1901575 1.201032e-04 2.139418e-02   260
GO:1990778 1.386091e-04 2.377006e-02    56
GO:0043085 1.452993e-04 2.377006e-02   160
GO:0044255 1.468243e-04 2.377006e-02   118
GO:0052547 1.474874e-04 2.377006e-02    64
GO:0006644 1.927152e-04 3.033696e-02    46
GO:0071705 2.101918e-04 3.201490e-02   273
GO:0006508 2.172251e-04 3.201490e-02   212
GO:0072659 2.175632e-04 3.201490e-02    50
GO:0055088 2.499678e-04 3.600068e-02    30
GO:0044281 2.935170e-04 4.077857e-02   204
GO:0051247 2.981310e-04 4.077857e-02   195
GO:0006753 3.012156e-04 4.077857e-02    53
GO:0033108 3.457867e-04 4.509012e-02    29
GO:0006629 3.463859e-04 4.509012e-02   173
GO:0030162 3.759149e-04 4.745525e-02   103
GO:0032268 3.785764e-04 4.745525e-02   302
GO:0009117 3.887757e-04 4.765326e-02    52
GO:0006793 4.041751e-04 4.765326e-02   307
GO:0008610 4.055773e-04 4.765326e-02    89
GO:0006886 4.152896e-04 4.765326e-02   138
GO:0006796 4.153557e-04 4.765326e-02   304
GO:0098876 4.280859e-04 4.809988e-02    22
GO:0006396 4.396619e-04 4.809988e-02   125
GO:0044093 4.405662e-04 4.809988e-02   207
GO:0090407 4.756379e-04 4.975830e-02    63
GO:0022904 4.778090e-04 4.975830e-02    28
GO:0042157 4.778090e-04 4.975830e-02    28
GO:0005739 6.429486e-13 5.117871e-10   298
GO:0031967 4.580767e-12 1.215430e-09   215
GO:0031975 4.580767e-12 1.215430e-09   215
GO:1902494 9.630289e-11 1.916428e-08   229
GO:0098798 5.208068e-10 8.291245e-08    83
GO:0005740 2.372151e-09 3.147054e-07   162
GO:0031966 1.677089e-08 1.907090e-06   150
GO:0005840 1.552868e-06 1.545104e-04    73
GO:0005743 3.425283e-06 2.834725e-04   108
GO:1990904 3.561212e-06 2.834725e-04   134
GO:0005730 4.548621e-06 3.291547e-04   140
GO:0098800 5.252726e-06 3.484308e-04    42
GO:0019866 7.817954e-06 4.673894e-04   112
GO:0044391 8.220417e-06 4.673894e-04    67
GO:0070062 2.330578e-05 1.236760e-03   244
GO:0005635 1.265685e-04 5.573592e-03    61
GO:0043230 1.330380e-04 5.573592e-03   244
GO:0065010 1.330380e-04 5.573592e-03   244
GO:1903561 1.330380e-04 5.573592e-03   244
GO:0140535 1.716557e-04 6.831017e-03   101
GO:0005783 1.802153e-04 6.831017e-03   268
GO:1990234 3.563484e-04 1.267772e-02   124
GO:0005759 3.822427e-04 1.267772e-02    90
GO:0048471 3.822427e-04 1.267772e-02    90
GO:1990204 7.008099e-04 2.231379e-02    27
GO:0098796 8.702866e-04 2.664416e-02   192
GO:0005746 1.025057e-03 3.022020e-02    31
GO:0005815 1.150636e-03 3.087819e-02    74
GO:0005856 1.187699e-03 3.087819e-02   233
GO:0005777 1.295821e-03 3.087819e-02    19
GO:0042579 1.295821e-03 3.087819e-02    19
GO:0000313 1.375887e-03 3.087819e-02    30
GO:0005761 1.375887e-03 3.087819e-02    30
GO:0098803 1.375887e-03 3.087819e-02    30
GO:0015934 1.396501e-03 3.087819e-02    39
[?25h

 > Analyzing the PROTECTED gene set 

[?25h[?25h           ONTOLOGY
GO:0007186       BP
GO:0007606       BP
GO:0050907       BP
GO:0009593       BP
GO:0007608       BP
GO:0050906       BP
GO:0050911       BP
GO:0007600       BP
GO:0051606       BP
GO:0050877       BP
GO:0003008       BP
GO:0007188       BP
GO:0007187       BP
GO:0007189       BP
GO:0048598       BP
GO:0048562       BP
GO:0050909       BP
GO:0007389       BP
GO:0007210       BP
GO:0098664       BP
GO:0007193       BP
GO:0007200       BP
GO:0050912       BP
GO:0050913       BP
GO:0007416       BP
GO:0006357       BP
GO:0051965       BP
GO:0001580       BP
GO:0001764       BP
GO:0003002       BP
GO:0009952       BP
GO:0051480       BP
GO:0006366       BP
GO:0051962       BP
GO:0007507       BP
GO:0048568       BP
GO:0007204       BP
GO:0051960       BP
GO:0042310       BP
GO:0051963       BP
GO:0072359       BP
GO:0009311       BP
GO:0003007       BP
GO:0034329       BP
GO:0006936       BP
GO:1903350       BP
GO:1903351       BP
GO:2000026       BP
GO:0003151       BP
GO:0009887       BP
GO:0009790       BP
GO:0007268       BP
GO:0098916       BP
GO:0071880       BP
GO:0090520       BP
GO:0050808       BP
GO:0099537       BP
GO:0019229       BP
GO:0048704       BP
GO:0007218       BP
GO:0009312       BP
GO:0005887       CC
GO:0031226       CC
GO:0005921       CC
GO:0005922       CC
GO:0000785       CC
GO:0099240       CC
GO:0099699       CC
GO:0043679       CC
GO:0004930       MF
GO:0004888       MF
GO:0038023       MF
GO:0060089       MF
GO:0004984       MF
GO:0001653       MF
GO:0008528       MF
GO:0003700       MF
GO:0000981       MF
GO:0000977       MF
GO:1990837       MF
GO:0043565       MF
GO:0003690       MF
GO:0000976       MF
GO:0001067       MF
GO:0005549       MF
GO:0140110       MF
GO:0000987       MF
GO:0000978       MF
GO:0003677       MF
GO:0008227       MF
GO:0001637       MF
GO:0004950       MF
GO:0008527       MF
GO:0030594       MF
GO:0016493       MF
GO:0033038       MF
GO:0019956       MF
GO:0019957       MF
GO:0042923       MF
GO:0015267       MF
                                                                                           Description
GO:0007186                                                G protein-coupled receptor signaling pathway
GO:0007606                                                     sensory perception of chemical stimulus
GO:0050907                               detection of chemical stimulus involved in sensory perception
GO:0009593                                                              detection of chemical stimulus
GO:0007608                                                                 sensory perception of smell
GO:0050906                                        detection of stimulus involved in sensory perception
GO:0050911                      detection of chemical stimulus involved in sensory perception of smell
GO:0007600                                                                          sensory perception
GO:0051606                                                                       detection of stimulus
GO:0050877                                                                      nervous system process
GO:0003008                                                                              system process
GO:0007188                   adenylate cyclase-modulating G protein-coupled receptor signaling pathway
GO:0007187 G protein-coupled receptor signaling pathway, coupled to cyclic nucleotide second messenger
GO:0007189                   adenylate cyclase-activating G protein-coupled receptor signaling pathway
GO:0048598                                                                     embryonic morphogenesis
GO:0048562                                                               embryonic organ morphogenesis
GO:0050909                                                                 sensory perception of taste
GO:0007389                                                               pattern specification process
GO:0007210                                                        serotonin receptor signaling pathway
GO:0098664                                      G protein-coupled serotonin receptor signaling pathway
GO:0007193                   adenylate cyclase-inhibiting G protein-coupled receptor signaling pathway
GO:0007200                     phospholipase C-activating G protein-coupled receptor signaling pathway
GO:0050912                      detection of chemical stimulus involved in sensory perception of taste
GO:0050913                                                          sensory perception of bitter taste
GO:0007416                                                                            synapse assembly
GO:0006357                                            regulation of transcription by RNA polymerase II
GO:0051965                                                     positive regulation of synapse assembly
GO:0001580               detection of chemical stimulus involved in sensory perception of bitter taste
GO:0001764                                                                            neuron migration
GO:0003002                                                                             regionalization
GO:0009952                                                    anterior/posterior pattern specification
GO:0051480                                           regulation of cytosolic calcium ion concentration
GO:0006366                                                          transcription by RNA polymerase II
GO:0051962                                           positive regulation of nervous system development
GO:0007507                                                                           heart development
GO:0048568                                                                 embryonic organ development
GO:0007204                                  positive regulation of cytosolic calcium ion concentration
GO:0051960                                                    regulation of nervous system development
GO:0042310                                                                            vasoconstriction
GO:0051963                                                              regulation of synapse assembly
GO:0072359                                                              circulatory system development
GO:0009311                                                           oligosaccharide metabolic process
GO:0003007                                                                         heart morphogenesis
GO:0034329                                                                      cell junction assembly
GO:0006936                                                                          muscle contraction
GO:1903350                                                                        response to dopamine
GO:1903351                                                               cellular response to dopamine
GO:2000026                                          regulation of multicellular organismal development
GO:0003151                                                                 outflow tract morphogenesis
GO:0009887                                                                  animal organ morphogenesis
GO:0009790                                                                          embryo development
GO:0007268                                                              chemical synaptic transmission
GO:0098916                                                        anterograde trans-synaptic signaling
GO:0071880                          adenylate cyclase-activating adrenergic receptor signaling pathway
GO:0090520                                                     sphingolipid mediated signaling pathway
GO:0050808                                                                        synapse organization
GO:0099537                                                                    trans-synaptic signaling
GO:0019229                                                              regulation of vasoconstriction
GO:0048704                                                     embryonic skeletal system morphogenesis
GO:0007218                                                              neuropeptide signaling pathway
GO:0009312                                                        oligosaccharide biosynthetic process
GO:0005887                                                       integral component of plasma membrane
GO:0031226                                                      intrinsic component of plasma membrane
GO:0005921                                                                                gap junction
GO:0005922                                                                            connexin complex
GO:0000785                                                                                   chromatin
GO:0099240                                                    intrinsic component of synaptic membrane
GO:0099699                                                     integral component of synaptic membrane
GO:0043679                                                                               axon terminus
GO:0004930                                                         G protein-coupled receptor activity
GO:0004888                                                   transmembrane signaling receptor activity
GO:0038023                                                                 signaling receptor activity
GO:0060089                                                               molecular transducer activity
GO:0004984                                                                 olfactory receptor activity
GO:0001653                                                                   peptide receptor activity
GO:0008528                                                 G protein-coupled peptide receptor activity
GO:0003700                                                   DNA-binding transcription factor activity
GO:0000981                       DNA-binding transcription factor activity, RNA polymerase II-specific
GO:0000977             RNA polymerase II transcription regulatory region sequence-specific DNA binding
GO:1990837                                               sequence-specific double-stranded DNA binding
GO:0043565                                                               sequence-specific DNA binding
GO:0003690                                                                 double-stranded DNA binding
GO:0000976                                                 transcription cis-regulatory region binding
GO:0001067                                        transcription regulatory region nucleic acid binding
GO:0005549                                                                             odorant binding
GO:0140110                                                            transcription regulator activity
GO:0000987                                         cis-regulatory region sequence-specific DNA binding
GO:0000978                       RNA polymerase II cis-regulatory region sequence-specific DNA binding
GO:0003677                                                                                 DNA binding
GO:0008227                                                   G protein-coupled amine receptor activity
GO:0001637                                         G protein-coupled chemoattractant receptor activity
GO:0004950                                                                 chemokine receptor activity
GO:0008527                                                                     taste receptor activity
GO:0030594                                                          neurotransmitter receptor activity
GO:0016493                                                             C-C chemokine receptor activity
GO:0033038                                                              bitter taste receptor activity
GO:0019956                                                                           chemokine binding
GO:0019957                                                                       C-C chemokine binding
GO:0042923                                                                        neuropeptide binding
GO:0015267                                                                            channel activity
                 pvalue     p.adjust Count
GO:0007186 1.902748e-57 1.000275e-53   240
GO:0007606 3.101574e-55 8.152487e-52   118
GO:0050907 3.113841e-52 5.456487e-49   110
GO:0009593 7.869513e-50 1.034251e-46   114
GO:0007608 4.272563e-49 4.492173e-46    99
GO:0050906 2.511081e-48 2.200125e-45   113
GO:0050911 1.926660e-46 1.446921e-43    94
GO:0007600 2.486316e-40 1.633820e-37   152
GO:0051606 1.205580e-39 7.041926e-37   120
GO:0050877 1.044980e-35 5.493459e-33   183
GO:0003008 1.105922e-28 5.285300e-26   230
GO:0007188 5.890483e-12 2.580522e-09    45
GO:0007187 7.921456e-12 3.203315e-09    28
GO:0007189 7.841907e-08 2.944636e-05    30
GO:0048598 8.643096e-08 3.029117e-05    62
GO:0048562 1.140276e-07 3.580607e-05    38
GO:0050909 1.157891e-07 3.580607e-05    19
GO:0007389 2.207064e-07 6.445854e-05    47
GO:0007210 8.159222e-07 2.144651e-04    14
GO:0098664 8.159222e-07 2.144651e-04    14
GO:0007193 9.110008e-07 2.206409e-04    17
GO:0007200 9.233592e-07 2.206409e-04    23
GO:0050912 1.620848e-06 3.550332e-04    15
GO:0050913 1.620848e-06 3.550332e-04    15
GO:0007416 2.312136e-06 4.861960e-04    25
GO:0006357 2.426815e-06 4.906833e-04   179
GO:0051965 3.396954e-06 6.613996e-04    17
GO:0001580 4.651156e-06 8.603768e-04    14
GO:0001764 4.746229e-06 8.603768e-04    19
GO:0003002 5.590917e-06 9.797151e-04    39
GO:0009952 8.779234e-06 1.488788e-03    27
GO:0051480 1.072475e-05 1.761875e-03    46
GO:0006366 1.645473e-05 2.621288e-03   181
GO:0051962 2.282073e-05 3.528488e-03    33
GO:0007507 3.459184e-05 5.149311e-03    46
GO:0048568 3.526254e-05 5.149311e-03    48
GO:0007204 4.709662e-05 6.515445e-03    42
GO:0051960 4.709662e-05 6.515445e-03    42
GO:0042310 6.033762e-05 8.133202e-03    14
GO:0051963 6.620359e-05 8.700807e-03    19
GO:0072359 9.205183e-05 1.180284e-02    76
GO:0009311 1.039017e-04 1.300503e-02    11
GO:0003007 1.068335e-04 1.306101e-02    26
GO:0034329 1.186550e-04 1.417658e-02    41
GO:0006936 1.447896e-04 1.649297e-02    33
GO:1903350 1.474548e-04 1.649297e-02    18
GO:1903351 1.474548e-04 1.649297e-02    18
GO:2000026 1.554970e-04 1.703016e-02   101
GO:0003151 1.637069e-04 1.756341e-02    14
GO:0009887 1.853025e-04 1.948270e-02    76
GO:0009790 2.100663e-04 2.165331e-02    78
GO:0007268 2.251163e-04 2.232899e-02    57
GO:0098916 2.251163e-04 2.232899e-02    57
GO:0071880 2.815700e-04 2.646349e-02     7
GO:0090520 2.815700e-04 2.646349e-02     7
GO:0050808 2.819014e-04 2.646349e-02    36
GO:0099537 2.936815e-04 2.708568e-02    57
GO:0019229 3.481670e-04 3.155714e-02    11
GO:0048704 3.899092e-04 3.474157e-02    14
GO:0007218 5.442028e-04 4.768124e-02    22
GO:0009312 5.693714e-04 4.906861e-02     8
GO:0005887 1.025102e-16 6.160862e-14   145
GO:0031226 3.459429e-15 1.039559e-12   149
GO:0005921 8.281161e-07 1.244244e-04    12
GO:0005922 8.281161e-07 1.244244e-04    12
GO:0000785 2.654872e-04 3.055274e-02    92
GO:0099240 4.025231e-04 3.055274e-02    14
GO:0099699 4.025231e-04 3.055274e-02    14
GO:0043679 4.066921e-04 3.055274e-02    13
GO:0004930 7.563515e-94 6.202083e-91   223
GO:0004888 2.627795e-75 1.077396e-72   237
GO:0038023 2.537990e-68 5.202880e-66   243
GO:0060089 2.537990e-68 5.202880e-66   243
GO:0004984 1.121588e-46 1.839405e-44    94
GO:0001653 2.354025e-17 2.757572e-15    46
GO:0008528 2.354025e-17 2.757572e-15    46
GO:0003700 3.091059e-16 3.168335e-14   143
GO:0000981 2.311735e-14 2.106248e-12   134
GO:0000977 1.740097e-13 1.426879e-11   134
GO:1990837 3.533008e-13 2.633697e-11   142
GO:0043565 5.567556e-13 3.804496e-11   147
GO:0003690 1.435498e-12 8.810630e-11   144
GO:0000976 1.611701e-12 8.810630e-11   138
GO:0001067 1.611701e-12 8.810630e-11   138
GO:0005549 1.725374e-11 8.842540e-10    21
GO:0140110 3.435263e-11 1.657009e-09   159
GO:0000987 1.586654e-09 7.155318e-08   113
GO:0000978 1.657939e-09 7.155318e-08   111
GO:0003677 2.630175e-09 1.027906e-07   173
GO:0008227 2.632442e-09 1.027906e-07    19
GO:0001637 6.491401e-09 2.314326e-07    16
GO:0004950 6.491401e-09 2.314326e-07    16
GO:0008527 2.120617e-08 7.245440e-07    15
GO:0030594 5.825490e-08 1.910761e-06    18
GO:0016493 6.922525e-08 2.102396e-06    14
GO:0033038 6.922525e-08 2.102396e-06    14
GO:0019956 2.258114e-07 6.172179e-06    13
GO:0019957 2.258114e-07 6.172179e-06    13
GO:0042923 2.258114e-07 6.172179e-06    13
GO:0015267 7.232003e-07 1.737211e-05    41
[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h
FALSE  TRUE 
 2115   328 
[?25h
FALSE  TRUE 
 2401    42 
[?25h

 >> Contingency tables 

[?25h     [,1] [,2]
TRUE  328   42
TRUE   46  120
[?25h[?25h     FALSE TRUE
[1,]  2401   42
[2,]   914  120
[?25h
	Fisher's Exact Test for Count Data

data:  rbind(table(protected$V1 %in% nonessential$gene), table(unprotected$V1 %in% nonessential$gene))
p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  5.189696 11.021757
sample estimates:
odds ratio 
  7.499359 

[?25h     FALSE TRUE
[1,]  2115  328
[2,]   988   46
[?25h
	Fisher's Exact Test for Count Data

data:  rbind(table(protected$V1 %in% common_essential$gene), table(unprotected$V1 %in% common_essential$gene))
p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.2136024 0.4140493
sample estimates:
odds ratio 
 0.3003023 

[?25h[?25h(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ l
ls: cannot open directory '.': No such file or directory
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ ll
ls: cannot open directory '.': No such file or directory
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ ll
ls: cannot open directory '.': No such file or directory
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ ls
ls: cannot open directory '.': No such file or directory
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ ..
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ l
[0m[01;32m00_geneLevelAnalysis.R[0m*  [01;32m02_statisticalAnalysis.R[0m*  [01;32m04_timingAnalysis.R[0m*      [01;34mdata[0m/     test.R
[01;32m01_binLevelAnalysis.R[0m*   [01;32m03_BRCA-PAAD_analysis.R[0m*   [01;32m05_protectionAnalysis.R[0m*  [01;34mresults[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ll
total 180K
drwxrwxrwx 1 laninst laninst 4,0K set  7 17:38 [0m[34;42m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  2 12:04 [01;34m..[0m/
-rwxr-xr-x 1 laninst laninst 8,7K set  6 18:36 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K ago 31 11:52 [01;34mdata[0m/
drwxrwxr-x 1 laninst laninst 4,0K ago 23 10:48 [01;34mresults[0m/
-rw-r--r-- 1 laninst laninst  43K set  7 17:38 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ll
total 180K
drwxrwxrwx 1 laninst laninst 4,0K set  7 17:38 [0m[34;42m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  2 12:04 [01;34m..[0m/
-rwxr-xr-x 1 laninst laninst 8,7K set  6 18:36 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K ago 31 11:52 [01;34mdata[0m/
drwxrwxr-x 1 laninst laninst 4,0K ago 23 10:48 [01;34mresults[0m/
-rw-r--r-- 1 laninst laninst  43K set  7 17:38 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ v[K..
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  2 12:04  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 17:38  [34;42mmutation_compensation[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ mkfit [K[K[K[Kfit[K[K[Kdi [Kt mutationCompensation

Command 'mkdit' not found, did you mean:

  command 'mkdir' from deb coreutils (8.30-3ubuntu2)
  command 'mkdic' from deb canna-utils (3.7p3-14)
  command 'mkdist' from deb libmodule-package-rdf-perl (0.014-1)

Try: apt install <deb name>

(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ mkdit mutationCompensation
[K(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ m

kdit mutationCompensation
[C[C[C[C[C
[K[A(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ mkdit mutationCompensation
[C[C[C
[K[A(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ mkdit mutationCompensation
[K(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ mkdit mutationCompensation[C[1P mutationCompensationr mutationCompensation
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ 
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  7 17:43  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 17:38  [34;42mmutation_compensation[0m/
drwxrwxr-x  1 laninst laninst 4,0K set  7 17:43  [01;34mmutationCompensation[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ cd mutation
bash: cd: mutation: No such file or directory
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ cd mutationCompensation/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 8,0K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:43 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:43 [01;34m..[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ gh repo clone fabio-alfieri/mutationCompensation
Cloning into 'mutationCompensation'...
remote: Enumerating objects: 6, done.[K
remote: Counting objects:  16% (1/6)[K
remote: Counting objects:  33% (2/6)[K
remote: Counting objects:  50% (3/6)[K
remote: Counting objects:  66% (4/6)[K
remote: Counting objects:  83% (5/6)[K
remote: Counting objects: 100% (6/6)[K
remote: Counting objects: 100% (6/6), done.[K
remote: Compressing objects:  33% (1/3)[K
remote: Compressing objects:  66% (2/3)[K
remote: Compressing objects: 100% (3/3)[K
remote: Compressing objects: 100% (3/3), done.[K
remote: Total 6 (delta 0), reused 0 (delta 0), pack-reused 0[K
Receiving objects:  16% (1/6)
Receiving objects:  33% (2/6)
Receiving objects:  50% (3/6)
Receiving objects:  66% (4/6)
Receiving objects:  83% (5/6)
Receiving objects: 100% (6/6)
Receiving objects: 100% (6/6), done.
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 12K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:45 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:43 [01;34m..[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:45 [01;34mmutationCompensation[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ..
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  7 17:43  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 17:38  [34;42mmutation_compensation[0m/
drwxrwxr-x  1 laninst laninst 4,0K set  7 17:45  [01;34mmutationCompensation[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ v
v: command not found
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ gh repo clone fabio-alfieri/mutationCompensation
fatal: destination path 'mutationCompensation' already exists and is not an empty directory.
exit status 128
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  7 17:43  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 17:38  [34;42mmutation_compensation[0m/
drwxrwxr-x  1 laninst laninst 4,0K set  7 17:45  [01;34mmutationCompensation[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ..
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~[00m$ ll
total 392K
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39 [0m[01;34m.[0m/
drwxr-xr-x 12 root    root    4,0K nov 16  2021 [01;34m..[0m/
-rw-r--r--  1 ieo5099 ieo5099 4,0K giu 14  2021 ._.
-rw-------  1 ieo5099 ieo5099  44K set  6 20:39 .bash_history
-rw-r--r--  1 ieo5099 ieo5099  220 apr 28  2021 .bash_logout
-rw-r--r--  1 ieo5099 ieo5099 4,7K nov  9  2021 .bashrc
drwxr-xr-x  2 ieo5099 ieo5099 4,0K lug 20  2021 [01;34mbin[0m/
drwx------ 15 ieo5099 ieo5099 4,0K giu  3 16:16 [01;34m.cache[0m/
drwxrwxr-x  2 ieo5099 ieo5099 4,0K giu 14  2021 [01;34m.conda[0m/
-rw-rw-r--  1 ieo5099 ieo5099   76 nov  8  2021 .condarc
drwx------ 16 ieo5099 ieo5099 4,0K giu 28 15:04 [01;34m.config[0m/
drwxrwxr-x  7 ieo5099 ieo5099 4,0K set  9  2021 [01;34m.cpan[0m/
drwxrwxr-x  8 ieo5099 ieo5099 4,0K ago 31 12:36 [01;34mDesktop2[0m/
drwxr-xr-x  1 laninst laninst 4,0K set  7 17:43 [01;34mDesktop_linux[0m/
-rw-r--r--  1 ieo5099 ieo5099 4,0K mag 17  2021 ._.DS_Store
-rw-r--r--  1 ieo5099 ieo5099 8,1K giu 14  2021 .DS_Store
drwxrwxr-x 15 ieo5099 ieo5099 4,0K set 10  2021 [01;34mensembl-vep[0m/
-rw-r--r--  1 ieo5099 ieo5099  151 giu  3 16:12 .gitconfig
drwx------  3 ieo5099 ieo5099 4,0K apr 28  2021 [01;34m.gnupg[0m/
drwxr-xr-x  4 ieo5099 ieo5099 4,0K apr 29  2021 [01;34m.local[0m/
drwxrwxr-x 20 ieo5099 ieo5099 4,0K giu 28 14:14 [01;34mminiconda3[0m/
drwxr-xr-x  1 laninst laninst 128K ago  1 16:06 [01;34mmountHD[0m/
drwxrwxr-x  5 ieo5099 ieo5099 4,0K feb 11  2022 [01;34mmskcc-vcf2maf-754d68a[0m/
drwx------  3 ieo5099 ieo5099 4,0K apr 28  2021 [01;34m.mysql[0m/
-rw-------  1 ieo5099 ieo5099   18 ago  4  2021 .mysql_history
-rw-r--r--  1 ieo5099 ieo5099  357 apr 28  2021 .pam_environment
-rw-r--r--  1 ieo5099 ieo5099 2,3K apr 22 19:59 PCAWG_producePlotsKaryotypes.R
drwxrwxr-x  5 ieo5099 ieo5099 4,0K set  9  2021 [01;34mperl5[0m/
-rw-r--r--  1 ieo5099 ieo5099  807 apr 28  2021 .profile
drwxr-xr-x  2 ieo5099 ieo5099 4,0K apr 28  2021 [01;34mPublic[0m/
-rw-------  1 ieo5099 ieo5099  44K mag 10  2021 .python_history
-rw-------  1 ieo5099 ieo5099    0 mag 10  2021 .python_history-53757.tmp
drwxrwxr-x  2 ieo5099 ieo5099 4,0K dic 23  2021 [01;34m.r[0m/
drwxr-xr-x  3 ieo5099 ieo5099 4,0K apr 28  2021 [01;34mR[0m/
-rw-r--r--  1 ieo5099 ieo5099  15K set  7 13:13 .Rhistory
drwxr-xr-x  6 ieo5099 ieo5099 4,0K lug 19 11:40 [01;34m.rstudio[0m/
drwx------  2 ieo5099 ieo5099 4,0K giu  3 15:52 [01;34m.ssh[0m/
drwxrwxr-x  7 ieo5099 ieo5099 4,0K set  9  2021 [01;34m.synapseCache[0m/
drwxrwxr-x  6 ieo5099 ieo5099 4,0K giu 16  2021 [01;34mvcftools[0m/
drwxrwxr-x  4 ieo5099 ieo5099 4,0K set  9  2021 [01;34m.vep[0m/
-rw-rw-r--  1 ieo5099 ieo5099  317 nov  8  2021 .wget-hsts
-rw-------  1 ieo5099 ieo5099  444 mar 28 16:42 .Xauthority
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~[00m$ cd Desktop_linux/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  7 17:43  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 17:38  [34;42mmutation_compensation[0m/
drwxrwxr-x  1 laninst laninst 4,0K set  7 17:45  [01;34mmutationCompensation[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ rm -r mutationCompensation/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll[K[Kv[Kgh repo clone fabio-alfieri/mutationCompensation
Cloning into 'mutationCompensation'...
remote: Enumerating objects: 6, done.[K
remote: Counting objects:  16% (1/6)[K
remote: Counting objects:  33% (2/6)[K
remote: Counting objects:  50% (3/6)[K
remote: Counting objects:  66% (4/6)[K
remote: Counting objects:  83% (5/6)[K
remote: Counting objects: 100% (6/6)[K
remote: Counting objects: 100% (6/6), done.[K
remote: Compressing objects:  33% (1/3)[K
remote: Compressing objects:  66% (2/3)[K
remote: Compressing objects: 100% (3/3)[K
remote: Compressing objects: 100% (3/3), done.[K
remote: Total 6 (delta 0), reused 0 (delta 0), pack-reused 0[K
Receiving objects:  16% (1/6)
Receiving objects:  33% (2/6)
Receiving objects:  50% (3/6)
Receiving objects:  66% (4/6)
Receiving objects:  83% (5/6)
Receiving objects: 100% (6/6)
Receiving objects: 100% (6/6), done.
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ LL
LL: command not found
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  7 17:45  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 17:38  [34;42mmutation_compensation[0m/
drwxrwxr-x  1 laninst laninst 4,0K set  7 17:45  [01;34mmutationCompensation[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ cd mutationCompensation/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 16K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:45 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:45 [01;34m..[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:45 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst  119 set  7 17:45 README.md
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ git pull [K[K[K[K[K[K[K[K[Kllù[K
total 136K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:46 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:45 [01;34m..[0m/
-rwxr-xr-x 1 laninst laninst 8,7K set  6 18:36 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:45 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst  119 set  7 17:45 README.md
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ f[Kgit pull
Already up-to-date.
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 136K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:46 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:45 [01;34m..[0m/
-rwxr-xr-x 1 laninst laninst 8,7K set  6 18:36 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst  119 set  7 17:45 README.md
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ git pull origin 
Already up-to-date.
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ rm *.T[KR
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 16K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:48 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:45 [01;34m..[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst  119 set  7 17:45 README.md
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ gitpu[K[K push
Username for 'https://github.com': ^C
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 16K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:48 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:45 [01;34m..[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst  119 set  7 17:45 README.md
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 16K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:48 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:45 [01;34m..[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst  119 set  7 17:45 README.md
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 136K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:48 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:45 [01;34m..[0m/
-rwxr-xr-x 1 laninst laninst 8,7K set  6 18:36 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst  119 set  7 17:45 README.md
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ .[Kll
total 136K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:48 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:45 [01;34m..[0m/
-rwxr-xr-x 1 laninst laninst 8,7K set  6 18:36 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst  119 set  7 17:45 README.md
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ git push
Username for 'https://github.com': LUCAspezia94!
Password for 'https://LUCAspezia94!@github.com': 
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ git push
Username for 'https://github.com': fabio-alfieri
Password for 'https://fabio-alfieri@github.com': 
remote: Support for password authentication was removed on August 13, 2021.
remote: Please see https://docs.github.com/en/get-started/getting-started-with-git/about-remote-repositories#cloning-with-https-urls for information on currently recommended modes of authentication.
fatal: Authentication failed for 'https://github.com/fabio-alfieri/mutationCompensation.git/'
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ git push
Username for 'https://github.com': fabio-alfieri
Password for 'https://fabio-alfieri@github.com': 
Everything up-to-date
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 136K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:48 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:45 [01;34m..[0m/
-rwxr-xr-x 1 laninst laninst 8,7K set  6 18:36 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst  119 set  7 17:45 README.md
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 136K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:48 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:45 [01;34m..[0m/
-rwxr-xr-x 1 laninst laninst 8,7K set  6 18:36 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst  119 set  7 17:45 README.md
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ..
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  7 17:45  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 17:38  [34;42mmutation_compensation[0m/
drwxrwxr-x  1 laninst laninst 4,0K set  7 17:48  [01;34mmutationCompensation[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ cd m[K [KmutationCompensation/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 136K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:48 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:45 [01;34m..[0m/
-rwxr-xr-x 1 laninst laninst 8,7K set  6 18:36 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst  119 set  7 17:45 README.md
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ git push
Everything up-to-date
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ git push origins[K my-new-branch
error: src refspec my-new-branch does not match any
[31merror: failed to push some refs to 'https://github.com/fabio-alfieri/mutationCompensation.git'
[m(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ..
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  7 17:45  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 17:38  [34;42mmutation_compensation[0m/
drwxrwxr-x  1 laninst laninst 4,0K set  7 17:48  [01;34mmutationCompensation[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ rm c[Kmt[KutationCompensation/
rm: cannot remove 'mutationCompensation/': Is a directory
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  7 17:45  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 17:38  [34;42mmutation_compensation[0m/
drwxrwxr-x  1 laninst laninst 4,0K set  7 17:48  [01;34mmutationCompensation[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ m[Krm -r mutaitonC[K[K[K[K[KtionCompensation/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  7 17:54  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 17:38  [34;42mmutation_compensation[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ git iti[Kit[K[K[K[K[K[K[K[Kmkdri[K[Kir myproject
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  7 17:55  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 17:38  [34;42mmutation_compensation[0m/
drwxrwxr-x  1 laninst laninst 4,0K set  7 17:55  [01;34mmyproject[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ cd myproject/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/myproject[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/myproject[00m$ ll
total 8,0K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:55 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:55 [01;34m..[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/myproject[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/myproject[00m$ git itin[K[K[Knit
[33mhint: Using 'master' as the name for the initial branch. This default branch name[m
[33mhint: is subject to change. To configure the initial branch name to use in all[m
[33mhint: of your new repositories, which will suppress this warning, call:[m
[33mhint: [m
[33mhint: 	git config --global init.defaultBranch <name>[m
[33mhint: [m
[33mhint: Names commonly chosen instead of 'master' are 'main', 'trunk' and[m
[33mhint: 'development'. The just-created branch can be renamed via this command:[m
[33mhint: [m
[33mhint: 	git branch -m <name>[m
Initialised empty Git repository in /home/ieo5099/Desktop_linux/myproject/.git/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/myproject[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/myproject[00m$ ll
total 12K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:55 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:55 [01;34m..[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:55 [01;34m.git[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/myproject[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/myproject[00m$ touch me[Knleson.txt
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/myproject[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/myproject[00m$ ll
total 12K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:56 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:55 [01;34m..[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:55 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst    0 set  7 17:56 mnleson.txt
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/myproject[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/myproject[00m$ git sa[Ktatus
On branch master

No commits yet

Untracked files:
  (use "git add <file>..." to include in what will be committed)
	[31mmnleson.txt[m

nothing added to commit but untracked files present (use "git add" to track)
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/myproject[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/myproject[00m$ git add
Nothing specified, nothing added.
[33mhint: Maybe you wanted to say 'git add .'?[m
[33mhint: Turn this message off by running[m
[33mhint: "git config advice.addEmptyPathspec false"[m
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/myproject[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/myproject[00m$ git add mnleson.txt 
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/myproject[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/myproject[00m$ ll
total 12K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:56 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:55 [01;34m..[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:56 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst    0 set  7 17:56 mnleson.txt
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/myproject[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/myproject[00m$ ..
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  7 17:55  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 17:38  [34;42mmutation_compensation[0m/
drwxrwxr-x  1 laninst laninst 4,0K set  7 17:56  [01;34mmyproject[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ cd[K[Kmkdir mutationCompensation
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  7 17:57  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 17:38  [34;42mmutation_compensation[0m/
drwxrwxr-x  1 laninst laninst 4,0K set  7 17:57  [01;34mmutationCompensation[0m/
drwxrwxr-x  1 laninst laninst 4,0K set  7 17:56  [01;34mmyproject[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux[00m$ cd mutaion[K[K[KtionCompensation/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 8,0K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:57 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:57 [01;34m..[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ mk[K[Kgit init
[33mhint: Using 'master' as the name for the initial branch. This default branch name[m
[33mhint: is subject to change. To configure the initial branch name to use in all[m
[33mhint: of your new repositories, which will suppress this warning, call:[m
[33mhint: [m
[33mhint: 	git config --global init.defaultBranch <name>[m
[33mhint: [m
[33mhint: Names commonly chosen instead of 'master' are 'main', 'trunk' and[m
[33mhint: 'development'. The just-created branch can be renamed via this command:[m
[33mhint: [m
[33mhint: 	git branch -m <name>[m
Initialised empty Git repository in /home/ieo5099/Desktop_linux/mutationCompensation/.git/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 12K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:57 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:57 [01;34m..[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:57 [01;34m.git[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 132K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:58 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:57 [01;34m..[0m/
-rwxr-xr-x 1 laninst laninst 8,7K set  6 18:36 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:57 [01;34m.git[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ git add 00_geneLevelAnalysis.R 00[K1_binLevelAnalysis.R 02_statisticalAnalysis.R 03_BRCA-PAAD_analysis.R 04_timingAnalysis.R 05_protectionAnalysis.R 0[K
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ l
[0m[01;32m00_geneLevelAnalysis.R[0m*  [01;32m02_statisticalAnalysis.R[0m*  [01;32m04_timingAnalysis.R[0m*
[01;32m01_binLevelAnalysis.R[0m*   [01;32m03_BRCA-PAAD_analysis.R[0m*   [01;32m05_protectionAnalysis.R[0m*
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 132K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:58 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:57 [01;34m..[0m/
-rwxr-xr-x 1 laninst laninst 8,7K set  6 18:36 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:58 [01;34m.git[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ git push
fatal: No configured push destination.
Either specify the URL from the command-line or configure a remote repository using

    git remote add <name> <url>

and then push using the remote name

    git push <name>

(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ git sa[Ktatus
On branch master

No commits yet

Changes to be committed:
  (use "git rm --cached <file>..." to unstage)
	[32mnew file:   00_geneLevelAnalysis.R[m
	[32mnew file:   01_binLevelAnalysis.R[m
	[32mnew file:   02_statisticalAnalysis.R[m
	[32mnew file:   03_BRCA-PAAD_analysis.R[m
	[32mnew file:   04_timingAnalysis.R[m
	[32mnew file:   05_protectionAnalysis.R[m

(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ git commit -m "add scirpts[K[K[K[K[Kripts"
[master (root-commit) 3f91a33] add scripts
 6 files changed, 2751 insertions(+)
 create mode 100755 00_geneLevelAnalysis.R
 create mode 100755 01_binLevelAnalysis.R
 create mode 100755 02_statisticalAnalysis.R
 create mode 100755 03_BRCA-PAAD_analysis.R
 create mode 100755 04_timingAnalysis.R
 create mode 100755 05_protectionAnalysis.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ git btamc[K[K[K[Kranch
[?1h=
* [32mmaster[m[m

[K[?1l>(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ git remote add oridin[K[K[Kgin https://github.com/fabio-alfieri/mutationCompensation.git
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutationCompensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutationCompensation[00m$ ll
total 132K
drwxrwxr-x 1 laninst laninst 4,0K set  7 17:58 [0m[01;34m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  7 17:57 [01;34m..[0m/
-rwxr-xr-x 1 laninst laninst 8,7K set  6 18:36 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  2o 31 12:36 [01;34mDesktop2[0m/
drwxr-xr-x  1 laninst laninst 4,0K set  8 11:18 [01;34mDesktop_linux[0m/
-rw-r--r--  1 ieo5099 ieo5099 4,0K mag 17  2021 ._.DS_Store
-rw-r--r--  1 ieo5099 ieo5099 8,1K giu 14  2021 .DS_Store
drwxrwxr-x 15 ieo5099 ieo5099 4,0K set 10  2021 [01;34mensembl-vep[0m/
-rw-r--r--  1 ieo5099 ieo5099  151 giu  3 16:12 .gitconfig
drwx------  3 ieo5099 ieo5099 4,0K apr 28  2021 [01;34m.gnupg[0m/
drwxr-xr-x  4 ieo5099 ieo5099 4,0K apr 29  2021 [01;34m.local[0m/
drwxrwxr-x 20 ieo5099 ieo5099 4,0K giu 28 14:14 [01;34mminiconda3[0m/
drwxr-xr-x  1 laninst laninst 128K ago  1 16:06 [01;34mmountHD[0m/
drwxrwxr-x  5 ieo5099 ieo5099 4,0K feb 11  2022 [01;34mmskcc-vcf2maf-754d68a[0m/
drwx------  3 ieo5099 ieo5099 4,0K apr 28  2021 [01;34m.mysql[0m/
-rw-------  1 ieo5099 ieo5099   18 ago  4  2021 .mysql_history
-rw-r--r--  1 ieo5099 ieo5099  357 apr 28  2021 .pam_environment
-rw-r--r--  1 ieo5099 ieo5099 2,3K apr 22 19:59 PCAWG_producePlotsKaryotypes.R
drwxrwxr-x  5 ieo5099 ieo5099 4,0K set  9  2021 [01;34mperl5[0m/
-rw-r--r--  1 ieo5099 ieo5099  807 apr 28  2021 .profile
drwxr-xr-x  2 ieo5099 ieo5099 4,0K apr 28  2021 [01;34mPublic[0m/
-rw-------  1 ieo5099 ieo5099  44K mag 10  2021 .python_history
-rw-------  1 ieo5099 ieo5099    0 mag 10  2021 .python_history-53757.tmp
drwxrwxr-x  2 ieo5099 ieo5099 4,0K dic 23  2021 [01;34m.r[0m/
drwxr-xr-x  3 ieo5099 ieo5099 4,0K apr 28  2021 [01;34mR[0m/
-rw-r--r--  1 ieo5099 ieo5099  15K set  7 13:13 .Rhistory
drwxr-xr-x  6 ieo5099 ieo5099 4,0K lug 19 11:40 [01;34m.rstudio[0m/
drwx------  2 ieo5099 ieo5099 4,0K giu  3 15:52 [01;34m.ssh[0m/
drwxrwxr-x  7 ieo5099 ieo5099 4,0K set  9  2021 [01;34m.synapseCache[0m/
drwxrwxr-x  6 ieo5099 ieo5099 4,0K giu 16  2021 [01;34mvcftools[0m/
drwxrwxr-x  4 ieo5099 ieo5099 4,0K set  9  2021 [01;34m.vep[0m/
-rw-rw-r--  1 ieo5099 ieo5099  317 nov  8  2021 .wget-hsts
-rw-------  1 ieo5099 ieo5099  444 mar 28 16:42 .Xauthority
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~[00m$ git co[Klone[K[K[K[K[K[K[K[K[Kcd Desktop_linux/u[Kmutation_compensation/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ll
total 304K
drwxrwxrwx 1 laninst laninst 4,0K set  8 12:11 [0m[34;42m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 11:18 [01;34m..[0m/
-rwxr-xr-x 1 laninst laninst 9,4K set  8 12:12 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K set  8 11:14 [01;34m06_proteinAggregation[0m/
drwxrwxr-x 1 laninst laninst 4,0K ago 31 11:52 [01;34mdata[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  8 12:12 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst 2,6K set  8 12:11 README.md
drwxrwxr-x 1 laninst laninst 4,0K ago 23 10:48 [01;34mresults[0m/
-rw-r--r-- 1 laninst laninst 155K set  8 12:13 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ cd [K[K[KRscir[K[Kritp[K[Kpt '[K00_geneLevelAnalysis.R 
[?25h[?25h[?25h[?25hError: please specify the analysis you want to perform!
Execution halted
[?25h(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ Rscript 00_geneLevelAnalysis.R 
[?25h[?25h[?25h[?25hUsage: 00_geneLevelAnalysis.R [options]


Options:
	-t CHARACTER, --tables=CHARACTER
		Options are: [y/n]

	-p CHARACTER, --plots=CHARACTER
		Options are: [y/n]

	-h, --help
		Show this help message and exit


Error: please specify the analysis you want to perform!
Execution halted
[?25h(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ Rscript 00_geneLevelAnalysis.R 
[?25h[?25h[?25h[?25h[?25h[?25h

This script 
	(1) estimates gene amplification frequency and mu score (takes several hours and cores);
	(2) produces gene-level correlations 

[ by default (1) is disables while (2) is running 
  set produce_tables = TRUE if you want (1) analysis ]

[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
^C
Warning messages:
1: In cor.test.default(x, y, method = "spearman") :
  Cannot compute exact p-value with ties
2: In cor.test.default(x, y, method = "spearman") :
  Cannot compute exact p-value with ties
Execution halted
[?25h(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ Rscript 00_geneLevelAnalysis.R 
[?25h[?25h[?25h[?25hUsage: 00_geneLevelAnalysis.R [options]


Options:
	-t TABLES, --tables=TABLES
		Options are: [y/n]

	-p PLOTS, --plots=PLOTS
		Options are: [y/n]

	-h, --help
		Show this help message and exit


Error: please specify the analysis you want to perform!
Execution halted
[?25h(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ Rscript 00_geneLevelAnalysis.R 
[?25h[?25h[?25h[?25hUsage: 00_geneLevelAnalysis.R [options]


Options:
	-t , --tables=
		Options are: [y/n]

	-p , --plots=
		Options are: [y/n]

	-h, --help
		Show this help message and exit


Error: please specify the analysis you want to perform!
Execution halted
[?25h(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ x[K[K(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ [K(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ [K(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ [K(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ Rscript 00_geneLevelAnalysis.R 
[?25h[?25h[?25h[?25h

 >> You chose default options: 
	(1) --tables 'n';
	(2) --statistics 'y' 



[?25h[?25h[?25h

 > This script 
	(1) estimates gene amplification frequency and mu score (takes several hours and cores);
	(2) produces gene-level correlations 


[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
     tumor_type type      condition                  corP               corS
cor          OV gene amplifications     0.115951533861878  0.107188642992291
cor1       BLCA gene amplifications    0.0478500663455604 0.0734670978402058
cor2       PCPG gene amplifications    0.0789169290576281 0.0847080630133289
cor3       PRAD gene amplifications   0.00232633349255485 0.0551826348909893
cor4       KIRC gene amplifications -0.000592896960146272 0.0141924417953718
cor5       MESO gene amplifications    0.0737354201110429 0.0910783136305054
cor6       TGCT gene amplifications     0.066125681849924 0.0469639110152404
                   p.corP               p.corS
cor   3.2703839736652e-29 3.52198325930251e-25
cor1 8.48397363425946e-09 8.88280718131685e-19
cor2  0.00352807324829162  0.00173395160329965
cor3    0.820693011455222 7.42733895672263e-08
cor4    0.954695371266016    0.173832991377281
cor5   0.0012908298728576 6.95678692828962e-05
cor6  0.00760903579235996   0.0581570270816269
Warning messages:
1: In cor.test.default(x, y, method = "spearman") :
  Cannot compute exact p-value with ties
2: In cor.test.default(x, y, method = "spearman") :
  Cannot compute exact p-value with ties
3: In cor.test.default(x, y, method = "spearman") :
  Cannot compute exact p-value with ties
4: In cor.test.default(x, y, method = "spearman") :
  Cannot compute exact p-value with ties
5: In cor.test.default(x, y, method = "spearman") :
  Cannot compute exact p-value with ties
6: In cor.test.default(x, y, method = "spearman") :
  Cannot compute exact p-value with ties
7: In cor.test.default(x, y, method = "spearman") :
  Cannot compute exact p-value with ties
[?25h

OUTPUT of the script: 
 	 (1) raw tables path: results/tables/00_geneLevelAnalysis/ 
[?25h	 (2) correlation plots path: results/plots/00_geneLevelAnalysis/ 
[?25h	 (3) gene-level correlations path: results/tables/02_produceStatistics/ 

[?25h[?25h(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ^C^C
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ 
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ^C
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ 
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ Rscript 00_geneLevelAnalysis.R 
[?25h[?25h[?25h[?25h

 >> You chose default options: 
	(1) --tables 'n';
	(2) --statistics 'y' 



[?25h[?25h[?25h

 > This script 
	(1) estimates gene amplification frequency and mu score (takes several hours and cores);
	(2) produces gene-level correlations 


[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h
 >> (2) analyis: produces gene-level correlations`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
     tumor_type type      condition                  corP               corS
cor          OV gene amplifications     0.115951533861878  0.107188642992291
cor1       BLCA gene amplifications    0.0478500663455604 0.0734670978402058
cor2       PCPG gene amplifications    0.0789169290576281 0.0847080630133289
cor3       PRAD gene amplifications   0.00232633349255485 0.0551826348909893
cor4       KIRC gene amplifications -0.000592896960146272 0.0141924417953718
cor5       MESO gene amplifications    0.0737354201110429 0.0910783136305054
cor6       TGCT gene amplifications     0.066125681849924 0.0469639110152404
                   p.corP               p.corS
cor   3.2703839736652e-29 3.52198325930251e-25
cor1 8.48397363425946e-09 8.88280718131685e-19
cor2  0.00352807324829162  0.00173395160329965
cor3    0.820693011455222 7.42733895672263e-08
cor4    0.954695371266016    0.173832991377281
cor5   0.0012908298728576 6.95678692828962e-05
cor6  0.00760903579235996   0.0581570270816269
[?25h

OUTPUT of the script: 
 	 (1) raw tables path: results/tables/00_geneLevelAnalysis/ 
[?25h	 (2) correlation plots path: results/plots/00_geneLevelAnalysis/ 
[?25h	 (3) gene-level correlations path: results/tables/02_produceStatistics/ 

[?25h[?25h(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ Rscript 00_geneLevelAnalysis.R 
[?25h[?25h[?25h[?25h

 >> You chose default options: 
	(1) --tables 'n';
	(2) --statistics 'y' 

[?25h[?25h[?25h

 > This script 
	(1) estimates gene amplification frequency and mu score (takes several hours and cores);
	(2) produces gene-level correlations 


[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h
 >> (2) analyis: produces gene-level correlations`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'


      tumor_type type      condition                  corP               corS
cor         LUAD gene amplifications    0.0815586605516272 0.0941689746648636
cor1        LUSC gene amplifications     0.049498551550676 0.0928482342201481
cor2        BRCA gene amplifications    0.0909297439881229   0.10483738418373
cor3        CESC gene amplifications    0.0766311273810531 0.0873952022819212
cor4        THCA gene amplifications    0.0544816941094508  0.113472394583287
cor5        HNSC gene amplifications    0.0105355714547944 0.0591454560461901
cor6        PAAD gene amplifications     0.083682450332083 0.0857227620103795
cor7    COADREAD gene amplifications     0.102286904312279  0.120218906852852
cor8      GBMLGG gene amplifications    0.0468363848480199  0.112508437970266
cor9          OV gene amplifications     0.115951533861878  0.107188642992291
cor10       BLCA gene amplifications    0.0478500663455604 0.0734670978402058
cor11       PCPG gene amplifications    0.0789169290576281 0.0847080630133289
cor12       PRAD gene amplifications   0.00232633349255485 0.0551826348909893
cor13       KIRC gene amplifications -0.000592896960146272 0.0141924417953718
cor14       MESO gene amplifications    0.0737354201110429 0.0910783136305054
cor15       TGCT gene amplifications     0.066125681849924 0.0469639110152404
cor16       KIRP gene amplifications     0.115554685062401 0.0419722490057031
cor17       SARC gene amplifications     0.056848269204535 0.0491379494782469
cor18       LIHC gene amplifications    0.0864731879960546 0.0915314827458046
cor19       ESCA gene amplifications    0.0559531789154683 0.0710943046356431
cor20       STAD gene amplifications    0.0213064024500785 0.0490855490358795
cor21        UCS gene amplifications     0.167423224804746  0.160782556180166
cor22       SKCM gene amplifications    0.0625940788365543 0.0469830071614574
                    p.corP               p.corS
cor   1.10924534199229e-22 1.00626460892131e-29
cor1  3.71200171586031e-09  1.6086342706167e-28
cor2   1.1694356258118e-26 6.23169626429976e-35
cor3   5.8320622485435e-19 3.34213332242429e-24
cor4   0.00126023300829389 1.65702407044605e-11
cor5      0.21682578414712 3.94611336516789e-12
cor6  4.48823521450506e-08 2.09045028603869e-08
cor7  5.36329538566715e-37 1.54626966647511e-50
cor8  3.83329920132751e-08 4.96884786203464e-40
cor9   3.2703839736652e-29 3.52198325930251e-25
cor10 8.48397363425946e-09 8.88280718131685e-19
cor11  0.00352807324829162  0.00173395160329965
cor12    0.820693011455222 7.42733895672263e-08
cor13    0.954695371266016    0.173832991377281
cor14   0.0012908298728576 6.95678692828962e-05
cor15  0.00760903579235996   0.0581570270816269
cor16 4.96291546250011e-28 7.02951091074589e-05
cor17 4.82119924692288e-07 1.36299359810838e-05
cor18 2.60464427908805e-20 1.42080301563749e-22
cor19 4.56345402614728e-08 3.64204987150906e-12
cor20  0.00932761568770448 2.06829202416857e-09
cor21 2.06759519892158e-29 3.17689833271859e-27
cor22 1.75893871625614e-11 4.52521012009249e-07
[?25h

OUTPUT of the script: 
 	 (1) raw tables path: results/tables/00_geneLevelAnalysis/ 
[?25h	 (2) correlation plots path: results/plots/00_geneLevelAnalysis/ 
[?25h	 (3) gene-level correlations path: results/tables/02_produceStatistics/ 

[?25h[?25h(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ Rscript 00_geneLevelAnalysis.R [KRscript 00_geneLevelAnalysis.R [1@-[1@-[1@v[1@a[1@n[1@i[1@l[1@l[1@a 00_geneLevelAnalysis.R  [A(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ [C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C
[?25h[?25h[?25h[?25h

 >> You chose default options: 
	(1) --tables 'n';
	(2) --statistics 'y' 

[?25h[?25h[?25h

 > This script 
	(1) estimates gene amplification frequency and mu score (takes several hours and cores);
	(2) produces gene-level correlations 


[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h
 >> (2) analyis: produces gene-level correlations`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
^C
Execution halted
[?25h(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ c[Kll
total 128K
drwxrwxrwx 1 laninst laninst 4,0K set  8 13:46 [0m[34;42m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 11:18 [01;34m..[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  8 11:14 [01;34m06_proteinAggregation[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 13:45 [01;34mcode[0m/
drwxrwxr-x 1 laninst laninst 4,0K ago 31 11:52 [01;34mdata[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  8 12:12 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst 2,6K set  8 12:11 README.md
drwxrwxr-x 1 laninst laninst 4,0K ago 23 10:48 [01;34mresults[0m/
-rw-r--r-- 1 laninst laninst  95K set  8 13:46 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ gir [K[Kty [K[K add code/[K[K[K[K[K[K[K[K[Kstatus
On branch main
Your branch is up-to-date with 'origin/main'.

Changes not staged for commit:
  (use "git add/rm <file>..." to update what will be committed)
  (use "git restore <file>..." to discard changes in working directory)
	[31mdeleted:    00_geneLevelAnalysis.R[m
	[31mdeleted:    01_binLevelAnalysis.R[m
	[31mdeleted:    02_statisticalAnalysis.R[m
	[31mdeleted:    03_BRCA-PAAD_analysis.R[m
	[31mdeleted:    04_timingAnalysis.R[m
	[31mdeleted:    05_protectionAnalysis.R[m

Untracked files:
  (use "git add <file>..." to include in what will be committed)
	[31mcode/[m
	[31mdata/[m
	[31mresults/[m
	[31mtest.R[m

no changes added to commit (use "git add" and/or "git commit -a")
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ git add coe[Kde/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ll
total 128K
drwxrwxrwx 1 laninst laninst 4,0K set  8 13:46 [0m[34;42m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 11:18 [01;34m..[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  8 11:14 [01;34m06_proteinAggregation[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 13:45 [01;34mcode[0m/
drwxrwxr-x 1 laninst laninst 4,0K ago 31 11:52 [01;34mdata[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  8 13:46 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst 2,6K set  8 12:11 README.md
drwxrwxr-x 1 laninst laninst 4,0K ago 23 10:48 [01;34mresults[0m/
-rw-r--r-- 1 laninst laninst  95K set  8 13:46 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ gfi[K[Kit stats
git: 'stats' is not a git command. See 'git --help'.

The most similar command is
	status
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ gia[Kt status
On branch main
Your branch is up-to-date with 'origin/main'.

Changes to be committed:
  (use "git restore --staged <file>..." to unstage)
	[32mnew file:   code/00_geneLevelAnalysis.R[m
	[32mnew file:   code/01_binLevelAnalysis.R[m
	[32mnew file:   code/02_statisticalAnalysis.R[m
	[32mnew file:   code/03_BRCA-PAAD_analysis.R[m
	[32mnew file:   code/04_timingAnalysis.R[m
	[32mnew file:   code/05_protectionAnalysis.R[m
	[32mnew file:   code/test.R[m

Changes not staged for commit:
  (use "git add/rm <file>..." to update what will be committed)
  (use "git restore <file>..." to discard changes in working directory)
	[31mdeleted:    00_geneLevelAnalysis.R[m
	[31mdeleted:    01_binLevelAnalysis.R[m
	[31mdeleted:    02_statisticalAnalysis.R[m
	[31mdeleted:    03_BRCA-PAAD_analysis.R[m
	[31mdeleted:    04_timingAnalysis.R[m
	[31mdeleted:    05_protectionAnalysis.R[m

Untracked files:
  (use "git add <file>..." to include in what will be committed)
	[31mdata/[m
	[31mresults/[m
	[31mtest.R[m

(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ git push
Username for 'https://github.com': fabio-alfieri
Password for 'https://fabio-alfieri@github.com': 
Everything up-to-date
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ git commint[K[Kt -m "ha[K[Kchange location to scrit[K[Kipts "
[main ee576aa] change location to scripts
 7 files changed, 4640 insertions(+)
 create mode 100755 code/00_geneLevelAnalysis.R
 create mode 100755 code/01_binLevelAnalysis.R
 create mode 100755 code/02_statisticalAnalysis.R
 create mode 100755 code/03_BRCA-PAAD_analysis.R
 create mode 100755 code/04_timingAnalysis.R
 create mode 100755 code/05_protectionAnalysis.R
 create mode 100644 code/test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ git puish[K[K[Ksh
Enumerating objects: 6, done.
Counting objects:  16% (1/6)Counting objects:  33% (2/6)Counting objects:  50% (3/6)Counting objects:  66% (4/6)Counting objects:  83% (5/6)Counting objects: 100% (6/6)Counting objects: 100% (6/6), done.
Delta compression using up to 48 threads
Compressing objects:  20% (1/5)Compressing objects:  40% (2/5)Compressing objects:  60% (3/5)Compressing objects:  80% (4/5)Compressing objects: 100% (5/5)Compressing objects: 100% (5/5), done.
Writing objects:  20% (1/5)Writing objects:  40% (2/5)Writing objects:  60% (3/5)Writing objects:  80% (4/5)Writing objects: 100% (5/5)Writing objects: 100% (5/5), 19.81 KiB | 1.80 MiB/s, done.
Total 5 (delta 1), reused 0 (delta 0), pack-reused 0
remote: Resolving deltas:   0% (0/1)[Kremote: Resolving deltas: 100% (1/1)[Kremote: Resolving deltas: 100% (1/1), completed with 1 local object.[K
To https://github.com/fabio-alfieri/mutation_compensation.git
   9e86d7e..ee576aa  main -> main
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ll
total 128K
drwxrwxrwx 1 laninst laninst 4,0K set  8 13:46 [0m[34;42m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 11:18 [01;34m..[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  8 11:14 [01;34m06_proteinAggregation[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 13:45 [01;34mcode[0m/
drwxrwxr-x 1 laninst laninst 4,0K ago 31 11:52 [01;34mdata[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  8 13:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst 2,6K set  8 12:11 README.md
drwxrwxr-x 1 laninst laninst 4,0K ago 23 10:48 [01;34mresults[0m/
-rw-r--r-- 1 laninst laninst  95K set  8 13:46 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ll
total 128K
drwxrwxrwx 1 laninst laninst 4,0K set  8 13:46 [0m[34;42m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 11:18 [01;34m..[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  8 11:14 [01;34m06_proteinAggregation[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 13:45 [01;34mcode[0m/
drwxrwxr-x 1 laninst laninst 4,0K ago 31 11:52 [01;34mdata[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  8 13:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst 2,6K set  8 12:11 README.md
drwxrwxr-x 1 laninst laninst 4,0K ago 23 10:48 [01;34mresults[0m/
-rw-r--r-- 1 laninst laninst  95K set  8 13:46 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ git push[K[K[K[K[K[K[K[Kgit [K[K[K[Kmv o[K06_proteinAggregation/
mv: missing destination file operand after '06_proteinAggregation/'
Try 'mv --help' for more information.
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ mv [K[K[Kll
total 128K
drwxrwxrwx 1 laninst laninst 4,0K set  8 13:46 [0m[34;42m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 11:18 [01;34m..[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  8 11:14 [01;34m06_proteinAggregation[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 13:45 [01;34mcode[0m/
drwxrwxr-x 1 laninst laninst 4,0K ago 31 11:52 [01;34mdata[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  8 13:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst 2,6K set  8 12:11 README.md
drwxrwxr-x 1 laninst laninst 4,0K ago 23 10:48 [01;34mresults[0m/
-rw-r--r-- 1 laninst laninst  95K set  8 13:46 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ mv 06_proteinAggregation/ code/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ l
[0m[01;34mcode[0m/  [01;34mdata[0m/  README.md  [01;34mresults[0m/  test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ll
total 124K
drwxrwxrwx 1 laninst laninst 4,0K set  8 13:48 [0m[34;42m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 11:18 [01;34m..[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 13:48 [01;34mcode[0m/
drwxrwxr-x 1 laninst laninst 4,0K ago 31 11:52 [01;34mdata[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  8 13:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst 2,6K set  8 12:11 README.md
drwxrwxr-x 1 laninst laninst 4,0K ago 23 10:48 [01;34mresults[0m/
-rw-r--r-- 1 laninst laninst  95K set  8 13:46 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ cd codeq
bash: cd: codeq: No such file or directory
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ l[Kcd code/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ ll
total 252K
drwxr-xr-x 1 laninst laninst 4,0K set  8 13:48 [0m[01;34m.[0m/
drwxrwxrwx 1 laninst laninst 4,0K set  8 13:48 [34;42m..[0m/
-rwxr-xr-x 1 laninst laninst 9,9K set  8 13:46 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K set  8 11:14 [01;34m06_proteinAggregation[0m/
-rw-r--r-- 1 laninst laninst 120K set  8 13:48 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ ..
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ll
total 124K
drwxrwxrwx 1 laninst laninst 4,0K set  8 13:48 [0m[34;42m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 11:18 [01;34m..[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 13:48 [01;34mcode[0m/
drwxrwxr-x 1 laninst laninst 4,0K ago 31 11:52 [01;34mdata[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  8 13:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst 2,6K set  8 12:11 README.md
drwxrwxr-x 1 laninst laninst 4,0K ago 23 10:48 [01;34mresults[0m/
-rw-r--r-- 1 laninst laninst  95K set  8 13:46 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ rm test.R 
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ll
total 28K
drwxrwxrwx 1 laninst laninst 4,0K set  8 13:48 [0m[34;42m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 11:18 [01;34m..[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  8 13:48 [01;34mcode[0m/
drwxrwxr-x 1 laninst laninst 4,0K ago 31 11:52 [01;34mdata[0m/
drwxrwxr-x 1 laninst laninst 4,0K set  8 13:47 [01;34m.git[0m/
-rw-rw-r-- 1 laninst laninst 2,6K set  8 12:11 README.md
drwxrwxr-x 1 laninst laninst 4,0K ago 23 10:48 [01;34mresults[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ cd code
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ ll
total 256K
drwxr-xr-x 1 laninst laninst 4,0K set  8 13:48 [0m[01;34m.[0m/
drwxrwxrwx 1 laninst laninst 4,0K set  8 13:48 [34;42m..[0m/
-rwxr-xr-x 1 laninst laninst 9,9K set  8 13:46 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K set  8 11:14 [01;34m06_proteinAggregation[0m/
-rw-r--r-- 1 laninst laninst 124K set  8 13:48 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ p[Kmore test.R 
#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(c("-a", "--analysis"), type="character", default=NULL, 
              help="Options are:
              \t - 00 or geneLevelAnalysis
              \t - 01 or binLevelAnalysis
              \t - 02 or produceStatistics
              \t - 03 or BRCAandPAADAnalysis
              \t - 04 or timingAnalysis
              \t - 05 or protectionAnalysis
              \t - all", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$analysis)) {
  print_help(opt_parser)
  stop("please specify the analysis you want to perform!", call.=FALSE)
}
if (!any(opt$analysis %in% c("00","geneLevelAnalysis",
                             "01","binLevelAnalysis",
                             "02","statisticalAnalysis",
                             "03","BRCAandPAADAnalysis",
                             "04","timingAnalysis",
                             "05","protectionAnalysis",
                             "all"
                             ))) {
  print_help(opt_parser)
  stop("typo in the analysis flag, plase see above for the available options!", call.=FALSE)
}

if (opt$analysis == "geneLevelAnalysis" | opt$analysis == "00") {
  system(paste0("Rscript --vanilla code/00_geneLevelAnalysis.R"))
}
[7m--More--(1%)[27m[K[7m--More--(1%)[27m[K
if (opt$analysis == "binLevelAnalysis" | opt$analysis == "01") {
  system(paste0("Rscript --vanilla code/01_binLevelAnalysis.R"))
}

if (opt$analysis == "protectionAnalysis" | opt$analysis == "05") {
  system(paste0("Rscript --vanilla code/05_protectionAnalysis.R"))
  stop()
}

# if (opt$analysis == "BRCAandPAADAnalysis" | opt$analysis == "03") {
#   system(paste0("Rscript --vanilla code/02_statisticalAnalysis.R"))
# }
# 
# if (opt$analysis == "statisticalAnalysis" | opt$analysis == "02") {
#   system(paste0("Rscript --vanilla code/02_statisticalAnalysis.R"))
# }
# 
# if (opt$analysis == "statisticalAnalysis" | opt$analysis == "02") {
#   system(paste0("Rscript --vanilla code/02_statisticalAnalysis.R"))
# }5h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[1mRows: [22m[34m1912[39m [1mColumns: [22m[34m1[39m
[36m──[39m [1mColumn specification[22m [36m─────────────────────────────────────────────────────────────────────────────────
───────────────────[39m
[1mDelimiter:[22m ","
[31mchr[39m (1): gene

[36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
[36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.
[?25h[1mRows: [22m[34m781[39m [1mColumns: [22m[34m1[39m
[36m──[39m [1mColumn specification[22m [36m─────────────────────────────────────────────────────────────────────────────────
───────────────────[39m
[1mDelimiter:[22m ","
[31mchr[39m (1): gene

[36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
[36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.
[?25h
[7m--More--(2%)[27mFALSE  TRUE 
 2115   328 
[?25h
FALSE  TRUE 
 2401    42 
[?25h

 >> Contingency tables 

[?25h     [,1] [,2]
TRUE  328   42
TRUE   46  120
[?25h[?25h     FALSE TRUE
[1,]  2401   42
[2,]   914  120
[?25h
	Fisher's Exact Test for Count Data

data:  rbind(table(protected$V1 %in% nonessential$gene), table(unprotected$V1 %in% nonessential$gene))
p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
  5.189696 11.021757
sample estimates:
odds ratio 
  7.499359 

[?25h     FALSE TRUE
[1,]  2115  328
[2,]   988   46
[?25h
	Fisher's Exact Test for Count Data

data:  rbind(table(protected$V1 %in% common_essential$gene), table(unprotected$V1 %in% common_essential$gene))
p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
[7m--More--(3%)[27m 0.2136024 0.4140493
sample estimates:
odds ratio 
 0.3003023 

[?25h[?25hError: 
Execution halted
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Worksta
tion[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ^C^C
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Worksta
tion[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ^C
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Worksta
tion[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ 
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Worksta
tion[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ^C
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Worksta
tion[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ^C
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Worksta
tion[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ 
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Worksta
tion[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ Rscript test.R -a 05[K[K[K[K[K[K[K[K[K[K[K[K
code/        data/        data.zip     results/     results.zip  Rplots.pdf   test.R       
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Worksta
tion[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ Rscript code/
00_geneLevelAnalysis.R    02_statisticalAnalysis.R  04_timingAnalysis.R       .Rapp.history
01_binLevelAnalysis.R     03_BRCA-PAAD_analysis.R   05_protectionAnalysis.R   
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Worksta
tion[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ Rscript code/05_protectionAnalysis.R 


 >> Loading libraries 

^C
Execution halted
[?25h(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Wo
rkstation[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ ..
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desk
[7m--More--(5%)[27mtop_linux[00m$ ll
total 2,9G
drwxr-xr-x  1 laninst laninst 4,0K set  2 12:04  [0m[01;34m.[0m/
drwxr-xr-x 25 ieo5099 ieo5099 4,0K set  6 20:39  [01;34m..[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m0_MAGISTRALE[0m/
drwxr-xr-x  1 laninst laninst 4,0K feb  4  2022  [01;34m8_TESI[0m/
-rw-rw-r--  1 laninst laninst  62K set  2 12:04  [01;35maggregation.png[0m
-rwxr-xr-x  1 laninst laninst 2,9G set 12  2021 [01;32m'LUCA PORTFOLIO.zip'[0m*
-rw-rw-r--  1 laninst laninst 4,7K lug 26  2021  maf2vcf.sh
drwxrwxr-x  1 laninst laninst 4,0K dic  3  2021  [01;34mmainProject_html[0m/
-rw-rw-r--  1 laninst laninst 261K ott  1  2021  [01;35mmeme.png[0m
drwxrwxr-x  1 laninst laninst 4,0K giu 15  2021  [01;34mMountWorkstation[0m/
drwxrwxrwx  1 laninst laninst 4,0K set  7 11:59  [34;42mmutation_compensation[0m/
-rw-rw-r--  1 laninst laninst  508 nov 12  2021  script.sh
-rw-rw-r--  1 laninst laninst 3,1K set  2 10:29  ubuntu.txt
drwxrwxr-x  1 laninst laninst 4,0K dic 13  2021  [01;34mxAngi[0m/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desk
top_linux[00m$ ll..Rscript code/05_protectionAnalysis.R 
Fatal error: cannot open file 'code/05_protectionAnalysis.R': No such file or directory
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desk
top_linux[00m$ cd mutation_compensation/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Worksta
tion[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ cd[K[Kll
total 6,6G
drwxrwxrwx 1 laninst laninst 4,0K set  7 11:59 [0m[34;42m.[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  2 12:04 [01;34m..[0m/
drwxr-xr-x 1 laninst laninst 4,0K set  6 17:59 [01;34mcode[0m/
drwxrwxr-x 1 laninst laninst 4,0K ago 31 11:52 [01;34mdata[0m/
-rw-rw-r-- 1 laninst laninst 5,4G set  1 16:38 [01;31mdata.zip[0m
drwxrwxr-x 1 laninst laninst 4,0K ago 23 10:48 [01;34mresults[0m/
-rw-rw-r-- 1 laninst laninst 1,3G set  1 16:25 [01;31mresults.zip[0m
-rw-rw-r-- 1 laninst laninst  24K set  7 17:32 Rplots.pdf
-rw-r--r-- 1 laninst laninst 6,1K set  7 17:34 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation[01;32mieo5099@laninst-HP-Z6-G4-Worksta
tion[00m:[01;34m~/Desktop_linux/mutation_compensation[00m$ cd code/
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Wo
rkstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ ll
[7m--More--(7%)[27mtotal 128K[K
drwxr-xr-x 1 laninst laninst 4,0K set  6 17:59 [0m[01;34m.[0m/
drwxrwxrwx 1 laninst laninst 4,0K set  7 11:59 [34;42m..[0m/
-rwxr-xr-x 1 laninst laninst 8,7K set  6 18:36 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:33 [01;32m05_protectionAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst    0 ago 30 15:03 [01;32m.Rapp.history[0m*
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Wo
rkstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ cd [K[K[KRscript --vai[Knilla 095[K[K0[K5_protectionAnalysis.R
 


 >> Loading libraries 

[?25h[?25h

 >> Calculate protection index 

[?25h[?25h		 LUAD 
Error in file(file, "rt") : cannot open the connection
Calls: suppressMessages ... withCallingHandlers -> read.csv -> read.table -> file
In addition: Warning message:
In file(file, "rt") :
  cannot open file 'results/tables/00_geneLevelAnalysis/LUAD_geneLevel.csv': No such file or directory
Execution halted
[?25h(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-
G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ Rscript --vanilla 05_protectionAnalysis.R 


 >> Loading libraries 

[?25h[?25h[?25h

 >> Calculate protection index 
[7m--More--(8%)[27m[K
[?25h[?25h		 LUAD 
		 LUSC 
		 BRCA 
		 CESC 
		 THCA 
		 HNSC 
		 PAAD 
		 GBMLGG 
		 COADREAD 
[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h

 >> Extacting common protected and unprotected genes across tumor types 

[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h

 >> Calculating essentiality scores for each gene of the two categories 

[?25h[?25h[?25h[?25h
	Welch Two Sample t-test

data:  mean_crispr_effect[mean_crispr_effect$genes %in% protected$V1, ]$mean and mean_crispr_effect[mean_crispr_effect$genes
 %in% unprotected$V1, ]$mean
t = -9.6362, df = 2647.7, p-value < 2.2e-16
alternative hypothesis: true difference in means is less than 0
95 percent confidence interval:
        -Inf -0.09674922
sample estimates:
  mean of x   mean of y 
-0.18701814 -0.07034669 

[?25h[?25h[?25h[?25h
	Wilcoxon rank sum test with continuity correction

data:  mean_crispr_effect[mean_crispr_effect$genes %in% protected$V1, ]$mean and mean_crispr_effect[mean_crispr_effect$genes
 %in% unprotected$V1, ]$mean
W = 828944, p-value < 2.2e-16
[7m--More--(9%)[27m[K(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ 
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ ll
total 268K
drwxr-xr-x 1 laninst laninst 4,0K set  8 13:48 [0m[01;34m.[0m/
drwxrwxrwx 1 laninst laninst 4,0K set  8 13:48 [34;42m..[0m/
-rwxr-xr-x 1 laninst laninst 9,9K set  8 13:46 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K set  8 11:14 [01;34m06_proteinAggregation[0m/
-rw-r--r-- 1 laninst laninst 136K set  8 13:49 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ ll
total 268K
drwxr-xr-x 1 laninst laninst 4,0K set  8 13:48 [0m[01;34m.[0m/
drwxrwxrwx 1 laninst laninst 4,0K set  8 13:48 [34;42m..[0m/
-rwxr-xr-x 1 laninst laninst 9,9K set  8 13:46 [01;32m00_geneLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  20K set  6 18:38 [01;32m01_binLevelAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  40K set  7 14:34 [01;32m02_statisticalAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  19K set  2 11:31 [01;32m03_BRCA-PAAD_analysis.R[0m*
-rwxr-xr-x 1 laninst laninst  13K set  1 17:53 [01;32m04_timingAnalysis.R[0m*
-rwxr-xr-x 1 laninst laninst  12K set  7 17:36 [01;32m05_protectionAnalysis.R[0m*
drwxrwxr-x 1 laninst laninst 4,0K set  8 11:14 [01;34m06_proteinAggregation[0m/
-rw-r--r-- 1 laninst laninst 136K set  8 13:49 test.R
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compensation/code[01;32mieo5099@laninst-HP-Z6-G4-Workstation[00m:[01;34m~/Desktop_linux/mutation_compensation/code[00m$ ..
(base) ]0;ieo5099@laninst-HP-Z6-G4-Workstation: ~/Desktop_linux/mutation_compens
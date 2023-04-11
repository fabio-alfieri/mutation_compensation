## Cancer genomes tolerate deleterious coding mutations through somatic copy-number amplifications of wild-type regions

For any questions regarding the code and/or data, please contact with Fabio Alfieri (fabio.alfieri@ieo.it).

[comment]: <> (Replace with the correct DOI)
[comment]:[![](https://img.shields.io/badge/doi-10.1101/2021.02.13.429885-rec.svg)](https://doi.org/10.1101/2021.03.13.429885)

### Abstract

Using 8,690 primary tumor samples from The Cancer Genome Atlas (TCGA) dataset, we performed an unbiased relationship analysis between somatic amplification and mutation frequency under different conditions and cancer types. We demonstrated that copy number amplifications frequently cover mutation-prone regions: this increases genomic tolerance towards the deleterious impact of mutations by creating extra safe copies of wild-type regions and hence protecting the genes. We showed that these potential compensatory events are highly influenced by gene and mutation properties, e.g., haploinsufficient genes showed a higher tendency to be protected from mutations via amplifications. Moreover, over-representation analysis of unprotected gene sets revealed unessential and cancer-unrelated function enrichment, as opposed to protected gene sets which contain cellular essential functions. A deeper understanding of these tumor genome-shaping compensatory events could reveal new cancer vulnerabilities and might guide novel therapeutic intervention.


### Reproduce analysis and figures

The available code was run and tested on Ubuntu (22.04.01 LTS), using R (4.2.1) and Python (3.9.12).

#### Clone the repository

```bash
git clone https://github.com/fabio-alfieri/mutation_compensation.git
```
#### Download supplementary data from Zenodo 

Download [here](https://doi.org/10.5281/zenodo.7079304) (Zenodo) the data.zip (required), put it in the cloned GitHub folder and unzip it:
```bash
cd path/to/GitHub/mutation_compensation/
wget -O data.zip https://zenodo.org/record/7079304/files/data.zip?download=1
unzip data.zip
```

To skip time- and computational-consuming steps, you can download the preprocessed outputs, the results.zip. In this case, put it in the cloned GitHub folder and unzip it:
```bash
cd path/to/GitHub/mutation_compensation/
wget -O results.zip https://zenodo.org/record/7079304/files/results.zip?download=1
unzip results.zip
```

#### Create conda environment

If you don't have conda installed locally, please [install Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html), and create the environment for the analysis:
```bash
cd path/to/GutHub/mutation_compensation/
conda env create -f conda/mutation_compensation.yml
```

#### Create folders

If you choose not to download the results.zip, please create the "results" folder and subfolders:
```bash
cd path/to/GitHub/mutation_compensation/
mkdir results
mkdir results/tables
mkdir results/plots
```

#### Run the analyses

In order to run the R scripts you should activate the conda environment and then launch the script:

```bash
conda activate mutation_compensation
cd path/to/GitHub/mutation_compensation/code/
Rscript 02_statisticalAnalysis.R --help
Rscript 02_statisticalAnalysis.R --tables y --statistics y
```

Scripts and detailed parameters:

| Rscript | Description | Parameters |
| --- | --- | --- |
| 00_geneLevelAnalysis.R | it produces gene level scores and correlations | `--tables` (y/[n]), produces gene-level score tables. If y, it may take several hours and it requires parallelization (see `cores` parameter within the Rscript and set according to your machine). Set n to skip this step, but only if you already downloaded the preprocessed data contained results.zip folder; `--statistics` ([y]/n), calculates gene level correlations |
| 01_binLevelAnalysis_TCGA.R | it produces bin level scores with multiple conditions using TCGA | no parameters |
| 01_binLevelAnalysis_PCAWG.R | it produces bin level scores with multiple conditions using PCAWG | no parameters |
| 02_statisticalAnalysis_TCGA.R | it produces correlations using different segmentation lengths and conditions | `--tables` ([y]/n), it may take several minutes depending on your machine; `--statistics` ([y]/n), produces statistics and plots |
| 02_statisticalAnalysis_PCAWG.R | it produces correlations using different segmentation lengths and conditions | no parameters |
| 03_BRCA-PAAD_analysis.R | it produces BRCA and PAAD analyses using OG and GO scores | no parameters |
| 04_timingAnalysis.R | it produces CNAqc and timing analyses on PCAWG dataset | `--runCNAqc` (y/[n]), if set to 'y' it may take several hours; `--test` (y/[n]), if set to 'y' it uses a small dataset (100,000 mutations) instead of the 5.5GB one |
| 04_bufferingAnalysis.R | it produces the buffering analyses on PCAWG dataset | no parameters |
| 05_protectionAnalysis.R | it etrieves the tables for protected and unprotected gene sets and perform Gene Ontology analyses |  `--runPermutations` (y/[n]), if set to 'y' it may take some time; it runs 10,000 permutations to caculate the distribution percentile of Pi scores |

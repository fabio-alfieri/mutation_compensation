## Cancer genomes tolerate deleterious coding mutations through somatic copy-number amplifications of wild-type regions

[comment]: <> (Replace with the correct DOI)
[comment]:[![](https://img.shields.io/badge/doi-10.1101/2021.02.13.429885-rec.svg)](https://doi.org/10.1101/2021.03.13.429885)

### Abstract

Using 8,690 primary tumor samples from The Cancer Genome Atlas (TCGA) dataset, we performed an unbiased relationship analysis between somatic amplification and mutation frequency under different conditions and cancer types. We demonstrated that copy number amplifications frequently cover mutation-prone regions: this increases genomic tolerance towards the deleterious impact of mutations by creating extra safe copies of wild-type regions and hence protecting the genes. We showed that these potential compensatory events are highly influenced by gene and mutation properties, e.g., haploinsufficient genes showed a higher tendency to be protected from mutations via amplifications. Moreover, over-representation analysis of unprotected gene sets revealed unessential and cancer-unrelated function enrichment, as opposed to protected gene sets which contain cellular essential functions. A deeper understanding of these tumor genome-shaping compensatory events could reveal new cancer vulnerabilities and might guide novel therapeutic intervention.

### Reproduce analysis and figures

#### (1) Clone this repository

Clone this repository:
```bash
git clone https://github.com/fabio-alfieri/mutation_compensation.git
```
#### (2) Dowload supplementary data from Zenodo 

Download from [here](https://doi.org/10.5281/zenodo.7065200) (Zenodo) the data.zip, put it in the cloned GitHub folder and unzip it:
```bash
cd path/to/GitHub/mutation_compensation/
wget -r https://zenodo.org/record/7065200/files/data.zip?download=1
unzip data.zip
```

#### (3) Dowload supplementary data from Zenodo 

If not, [install Conda](https://docs.conda.io/projects/conda/en/latest/commands/install.html), and create the environment for the analysis:
```bash
cd path/to/GutHub/mutation_compensation/
conda env create -f conda/mutation_compensation.yml
```

#### (4) Run create the folders

Create results folders and subfolders
```bash
cd path/to/GitHub/mutation_compensation/
mkdir results
mkdir results/tables
mkdir results/plots
```

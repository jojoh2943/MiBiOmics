# MiBiOmics
MiBiOmics is a shiny-based web application to perform correlation network analysis and multi-omics analysis on omics datasets

## Installation

MiBiOmics will soon be available at https://mibiomics.univ-nantes.fr/ but can also be runned locally in RStudio. To run MiBiOmics locally you need to install the lastest version of shiny, **R (version 3.5)** and run the following commands:

```r
library(shiny)
runGitHub("MiBiOmics", "jojoh2943")
```
All the necessary packages will be installed automatically for the usage of MiBiOmics. If you don't have the right R and shiny version you can download the conda environment from the file *MiBiOmics.yml* and run Rstudio from this environment. If your using the conda environment, you will have to run the following command in conda:

```bash
conda env create -f MiBiOmics.yml
conda activate MiBiOmics
R -e 'shiny::runGitHub("MiBiOmics", "jojoh2943")'
```

The first launch will take some time (about one hour) due to the installation of all the necessary packages.

## Data Upload & pre-treatment:

In the first section you can upload one, two or even three omics tables:

|        | Gene1  | Gene2 | Gene3 |
|        | ------ | ----- | ----- |
|Sample 1|   1    |   0   |  512  |
|Sample 2|  106   |   0   |   26  |
|Sample 3|   0    |   3   |   0   |
|Sample 4|   19   |   0   |   6   |

Omics table may contain Gene names, OTUs/ASVs or even metabolites in columns. Rownames must be unique. 
An annotation table describing your samples is also necessary:


|        |  Site  |  ADN  | group |
|        | ------ | ----- | ----- |
|Sample 1|   1    |  0.63 |   A   |
|Sample 2|   2    |  1.12 |   B   |
|Sample 3|   2    |  0.45 |   A   |
|Sample 4|   5    |  1.04 |   A   |

If an Omics table containing ASVs or OTUs is uploaded, you can also add a taxonomic table describing the phylogenetic information of each OTUs/ASVs.

|        |Kingdom |Phylum | Class | Order |Family |Species|
|        | ------ | ----- | ----- | ----- | ----- | ----- |
| OTU 1  |  ...   |  ...  |  ...  |  ...  |  ...  |  ...  |
| OTU 2  |  ...   |  ...  |  ...  |  ...  |  ...  |  ...  |
| OTU 3  |  ...   |  ...  |  ...  |  ...  |  ...  |  ...  |
| OTU 4  |  ...   |  ...  |  ...  |  ...  |  ...  |  ...  |


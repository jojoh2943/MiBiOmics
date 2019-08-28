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

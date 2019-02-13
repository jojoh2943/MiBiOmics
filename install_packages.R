
bioconductor_packages <- c("GO.db", "preprocessCore", "impute",
                           "sva", "metagenomeSeq", "omicade4",
                           "threejs")


install_packages <- c("shiny", "matrixStats", "Hmisc",
                      "splines", "foreach", "doParallel",
                      "fastcluster", "dynamicTreeCut", "survival",
                      "WGCNA", "flexdashboard", "plotly",
                      "shinythemes", "ggplot2",
                      "ggrepel", "ggdendro", "ggthemes",
                      "dendextend", "pheatmap", "reshape2",
                      "gridExtra", "RColorBrewer", "vegan",
                      "DT", "shinyjs", "cowplot",
                      "visNetwork", "GGally", "sna",
                      "pls", "BiocManager", "iheatmapr",
                      "rmarkdown", "plotly", "ade4",
                      "igraph", "network", "plotly",
                      "compositions")


biocManager_packages <- c("mixOmics")

if (length(setdiff(bioconductor_packages, rownames(installed.packages()))) > 0) {
  source("https://bioconductor.org/biocLite.R")
  biocLite(setdiff(bioconductor_packages, rownames(installed.packages())), suppressUpdates = TRUE)
}

if (length(setdiff(install_packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(install_packages, rownames(installed.packages())))
}


if (length(setdiff(biocManager_packages, rownames(installed.packages()))) > 0) {
  source("http://bioconductor.org/biocLite.R")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    +     install.packages("BiocManager")
  BiocManager::install("mixOmics", version = "3.8")
}


#define MAX_NUM_DLLS 10000
library("flexdashboard")
library("metagenomeSeq")
library("shiny")
library("V8")
library("shinythemes")
library("WGCNA")
library("ggplot2")
library("ggdendro")
library("ggrepel")
library('dendextend')
library("ggthemes")
library("pheatmap")
library("reshape2")
library("gridExtra")
library("RColorBrewer")
library("dynamicTreeCut")
library("omicade4")
library("ade4")
library("data.table")
library("vegan")
library("tools")
library("compositions")
library("DT")
library("plotly")
library("iheatmapr")
library("threejs")
library("shinyjs")
library("cowplot")
library("igraph")
library("visNetwork")
library("GGally")
library("network")
library("rmarkdown")
library("pls")
library("sva") #Batch correction
library("mixOmics") #CLR transformation
library("plotly")

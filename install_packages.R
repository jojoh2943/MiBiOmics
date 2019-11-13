
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
                      "compositions", "leaflet", "grid", "webshot", "igraph", "psych", "htmltools")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (length(setdiff(c(bioconductor_packages, install_packages), rownames(installed.packages()))) > 0) {
  BiocManager::install(setdiff(c(bioconductor_packages, install_packages), rownames(installed.packages())), updates = FALSE)
}

i <- installed.packages()

# prevent error writeImpl --> dur to version 0.4.0 of htmltools loaded as a dependency of leaflet
if ("htmltools" %in% rownames(installed.packages())){
  i <- installed.packages()
  if (i["htmltools", "Version"] != "0.3.6"){
    remove.packages("htmltools")
    packageurl <- "https://cran.r-project.org/src/contrib/Archive/htmltools/htmltools_0.3.6.tar.gz"
    install.packages(packageurl, repos=NULL, type="source")
  }
}

#define MAX_NUM_DLLS 10000
library("flexdashboard")
library("metagenomeSeq")
library("shiny")
#library("V8")
library("shinythemes")
library("WGCNA")
library("ggplot2")
library("ggdendro")
library("ggrepel")
library("dendextend")
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
library("plotly")
library("leaflet")
library("grid")
library("webshot") # save iheatmap
library("igraph") # Keystone index
library("psych") # factor analysis




#############################################
#### SERVER FUNCTIONS : DATA EXPLORATION ####
#############################################
# In this script you will find the functions used in the second tabulation of MiBiOmics called data exploration
# The goal of these functions is to explore the general data using dendrogrammes, PCoA, PCA and when working
# with microbiota datasets, relative abundance plots.


#### REACTIVE OBJECTS ####

#### DENDROGRAMME
dendrogramme <- reactive({
  sampleTree = hclust(dist(exprDat_2()), method = input$selectMethod)
  sampleTree
})

# Dendrogramme of the second dataset
dendrogramme_sec <- reactive({ 
  sampleTree = hclust(dist(exprDatSec_2()), method = input$selectMethod_sec)
  sampleTree
})



info_dendrogramme <- reactive({
  ddata_dendro <- dendro_data(dendrogramme())
  labs <- label(ddata_dendro)
  Samples = labs$label
  conditionRow = match(Samples, rownames(sampleAnnot_2()))
  dataCondition = sampleAnnot_2()[conditionRow, ]
  labs <- data.frame(labs, dataCondition)
  colnames(labs) <- c("x", "y", "label", colnames(dataCondition))
  
  labs
})


# Info Dendrogramme for the second dataset
info_dendrogramme_sec <- reactive({
  ddata_dendro <- dendro_data(dendrogramme_sec())
  labs <- label(ddata_dendro)
  Samples = labs$label
  conditionRow = match(Samples, rownames(sampleAnnot_2()))
  dataCondition = sampleAnnot_2()[conditionRow, ]
  labs <- data.frame(labs, dataCondition)
  colnames(labs) <- c("x", "y", "label", colnames(dataCondition))
  
  labs
})

#### RELATIVE ABUNDANCE PLOTS

exprDat_Relative_abundance <- reactive({
  validate(
    need(input$LoadExample == "Yes" || input$LoadExample2 == "Yes" || (input$TaxonFile == TRUE|| input$TaxonFile1 == TRUE), "This plot can only be realised with OTUs")
  )
  
  exprDat <- exprDat_2()
  
  if (input$SelectVariable_rel_ab != "unordered"){
    exprDat <- exprDat[order(sampleAnnot_2()[,input$SelectVariable_rel_ab]),]
  }
  exprDat_rel_ab <- as.data.frame(t(apply(exprDat, 1, function(x) x / sum(x))))
  exprDat_rel_ab <- setDT(as.data.frame(t(exprDat_rel_ab)), keep.rownames = TRUE)
  exprDat_rel_ab <- merge(exprDat_rel_ab, taxTable1(), by = "rn")
  exprDat_rel_ab$rn <- as.character(exprDat_rel_ab$rn)
  exprDat_rel_ab$rn <- factor(exprDat_rel_ab$rn, levels = unique(exprDat_rel_ab$rn))
  df_long <- melt(exprDat_rel_ab)
  df_long
})

exprDat_Relative_abundance_sec <- reactive({
  validate(
    need((input$OmicTable == "OTUs"), "This plot can only be realised with OTUs")
  )
  exprDat <- exprDatSec_3()
  
  if (input$SelectVariable_rel_ab_sec != "unordered"){
    exprDat <- exprDat[order(sampleAnnot_2()[,input$SelectVariable_rel_ab_sec]),]
  }
  exprDat_rel_ab <- as.data.frame(t(apply(exprDat, 1, function(x) x / sum(x))))
  exprDat_rel_ab <- setDT(as.data.frame(t(exprDat_rel_ab)), keep.rownames = TRUE)
  exprDat_rel_ab <- merge(exprDat_rel_ab, taxTable1(), by = "rn")
  exprDat_rel_ab$rn <- as.character(exprDat_rel_ab$rn)
  exprDat_rel_ab$rn <- factor(exprDat_rel_ab$rn, levels = unique(exprDat_rel_ab$rn))
  df_long <- melt(exprDat_rel_ab)
  df_long
})

#### ORDINATIONS

#### PCoA
dudi_PCoA <- reactive({
  if (input$selectOrdination == "PCoA"){
    OTU_dist <- vegdist(exprDat_2(), method = input$selectDist)
    dudi_Pcoa_OTU <- dudi.pco(OTU_dist, scannf = FALSE, nf = 2)
  }else{
    dudi_Pcoa_OTU <- dudi.pca(exprDat_2(), scannf = FALSE, nf = 2)
  }
  
  dudi_Pcoa_OTU
})

dudi_PCoA_sec <- reactive({
  if (input$selectOrdinationSec == "PCoA"){
    dudi_Pcoa_OTU <- dudi.pca(exprDatSec_3(), scannf = FALSE, nf = 2)
  }else{
    OTU_dist <- vegdist(exprDatSec_3(), method = input$selectDist_sec)
    dudi_Pcoa_OTU <- dudi.pco(OTU_dist, scannf = FALSE, nf = 2)
  }
  
  dudi_Pcoa_OTU
})

PCoA <- reactive({
  sample_info <- data.frame(dudi_PCoA()$li[,1:2], sampleAnnot_2())
  sample_info
})

PCoA_sec <- reactive({
  sample_info <- data.frame(dudi_PCoA_sec()$li[,1:2], sampleAnnot_sec())
  sample_info
})


#### OBSERVE EVENTS ####

observeEvent(input$showMethod, {
  showModal(modalDialog(
    title = "How to choose the right method ?",
    HTML("- Ward's minimum variance method aims at finding compact, spherical clusters. <br> <br>
         - The complete linkage method finds similar clusters. <br> <br>
         - The single linkage method (which is closely related to the minimal spanning tree) adopts a âfriends of friendsâ clustering strategy. <br> <br>
         - The other methods can be regarded as aiming for clusters with characteristics somewhere between the single and complete link methods."),
    easyClose = TRUE,
    size = "l"
    ))
})

observeEvent(input$showMethod_sec, {
  showModal(modalDialog(
    title = "How to choose the right method ?",
    HTML("- Ward's minimum variance method aims at finding compact, spherical clusters. <br> <br>
         - The complete linkage method finds similar clusters. <br> <br>
         - The single linkage method (which is closely related to the minimal spanning tree) adopts a âfriends of friendsâ clustering strategy. <br> <br>
         - The other methods can be regarded as aiming for clusters with characteristics somewhere between the single and complete link methods."),
    easyClose = TRUE,
    size = "l"
    ))
})

observeEvent(input$showMethod_P6, {
  showModal(modalDialog(
    title = "How to choose the right method ?",
    HTML("- Ward's minimum variance method aims at finding compact, spherical clusters. <br> <br>
         - The complete linkage method finds similar clusters. <br> <br>
         - The single linkage method (which is closely related to the minimal spanning tree) adopts a âfriends of friendsâ clustering strategy. <br> <br>
         - The other methods can be regarded as aiming for clusters with characteristics somewhere between the single and complete link methods."),
    easyClose = TRUE,
    size = "l"
    ))
})



observeEvent(input$showDist, {
  showModal(modalDialog(
    title = "How to choose the right distance ?",
    HTML("- Binomial index is derived from Binomial deviance under null hypothesis that the two compared communities are equal. It should be able to handle variable sample sizes. The index does not have a fixed upper limit, but can vary among sites with no shared species.<br> <br> - Cao index or CYd index (Cao et al. 1997) was suggested as a minimally biased index for high beta diversity and variable sampling intensity. Cao index does not have a fixed upper limit, but can vary among sites with no shared species. The index is intended for count (integer) data, and it is undefined for zero abundances; these are replaced with arbitrary value 0.1 <br> <br> - Euclidean and Manhattan dissimilarities are not good in gradient separation without proper standardization <br> <br> - BrayâCurtis and Jaccard indices are rank-order similar, and some other indices become identical or rank-order similar after some standardizations. Jaccard index is metric, and probably should be preferred instead of the default Bray-Curtis which is semimetric"),
    easyClose = TRUE,
    size = "l"
  ))
})

observeEvent(input$showDist_sec, {
  showModal(modalDialog(
    title = "How to choose the right distance ?",
    HTML("- Binomial index is derived from Binomial deviance under null hypothesis that the two compared communities are equal. It should be able to handle variable sample sizes. The index does not have a fixed upper limit, but can vary among sites with no shared species.<br> <br> - Cao index or CYd index (Cao et al. 1997) was suggested as a minimally biased index for high beta diversity and variable sampling intensity. Cao index does not have a fixed upper limit, but can vary among sites with no shared species. The index is intended for count (integer) data, and it is undefined for zero abundances; these are replaced with arbitrary value 0.1 <br> <br> - Euclidean and Manhattan dissimilarities are not good in gradient separation without proper standardization <br> <br> - BrayâCurtis and Jaccard indices are rank-order similar, and some other indices become identical or rank-order similar after some standardizations. Jaccard index is metric, and probably should be preferred instead of the default Bray-Curtis which is semimetric"),
    easyClose = TRUE,
    size = "l"
  ))
})

observeEvent(input$showDist_P6, {
  showModal(modalDialog(
    title = "How to choose the right distance ?",
    HTML("- Binomial index is derived from Binomial deviance under null hypothesis that the two compared communities are equal. It should be able to handle variable sample sizes. The index does not have a fixed upper limit, but can vary among sites with no shared species.<br> <br> - Cao index or CYd index (Cao et al. 1997) was suggested as a minimally biased index for high beta diversity and variable sampling intensity. Cao index does not have a fixed upper limit, but can vary among sites with no shared species. The index is intended for count (integer) data, and it is undefined for zero abundances; these are replaced with arbitrary value 0.1 <br> <br> - Euclidean and Manhattan dissimilarities are not good in gradient separation without proper standardization <br> <br> - BrayâCurtis and Jaccard indices are rank-order similar, and some other indices become identical or rank-order similar after some standardizations. Jaccard index is metric, and probably should be preferred instead of the default Bray-Curtis which is semimetric"),
    easyClose = TRUE,
    size = "l"
  ))
})

#### INTERFACE VARIABLES ####

output$SelectVariable <-  renderUI({
  selectInput("selectVariable", 
              label = "Choose a variable: ", 
              choices = colnames(sampleAnnot_2()), 
              selected = colnames(sampleAnnot_2())[2])
})

output$SelectVariable3 <-  renderUI({
  selectInput("selectVariable3", 
              label = "Choose a variable: ", 
              choices = colnames(sampleAnnot_2()), 
              selected = colnames(sampleAnnot_2())[2])
})

output$SelectTaxo <- renderUI({
  selectInput("selectTaxo",
              label = "Select a taxonomic level for the relative abundance plot",
              choices = colnames(taxTable_present()),
              selected = colnames(taxTable_present())[3])
})

output$SelectVariable_sec <-  renderUI({
  selectInput("selectVariable_sec", 
              label = "Choose a variable: ", 
              choices = colnames(sampleAnnot_2()), 
              selected = colnames(sampleAnnot_2())[2])
})

output$SelectVariable3_sec <-  renderUI({
  selectInput("selectVariable3_sec", 
              label = "Choose a variable: ", 
              choices = colnames(sampleAnnot_2()), 
              selected = colnames(sampleAnnot_2())[2])
})

output$SelectTaxo_sec <- renderUI({
  selectInput("selectTaxo_sec",
              label = "Select a taxonomic level for the relative abundance plot",
              choices = colnames(taxTable_present()),
              selected = colnames(taxTable_present())[3])
})

output$SelectVariable_rel_ab <- renderUI({
  selectInput("SelectVariable_rel_ab", 
              label = "Choose a variable: ", 
              choices = c("unordered", colnames(sampleAnnot_2())), 
              selected = "unordered")
})

output$SelectVariable_rel_ab_sec <- renderUI({
  selectInput("SelectVariable_rel_ab_sec", 
              label = "Choose a variable: ", 
              choices = c("unordered", colnames(sampleAnnot_2())), 
              selected = "unordered")
})




#### PLOT OUTPUTS ####

# General Relative Abundance
output$Relative_Abundance <- renderPlot({
  ggplot(exprDat_Relative_abundance(), aes_string(x = "variable", y = "value", fill = input$selectTaxo)) +
    geom_bar(stat = "identity") +
    labs(x = "Samples", 
         y = "Relative Abundance", 
         title = paste("Relative abundance at the ", input$selectTaxo, " Taxonomic level", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
})

# Dendrogramme
output$Dendrogramme <- renderPlot({
  if (!is.numeric(sampleAnnot_2()[,input$SelectVariable])){
    ggdendrogram(dendrogramme(), label = FALSE) +
      geom_text(data = info_dendrogramme(),
                aes_string(label = colnames(info_dendrogramme()[3]), x = colnames(info_dendrogramme()[1]), y = 0, colour = input$selectVariable),
                angle = 75,
                vjust = 0,
                nudge_y = -0.1,
                size = 3) +
      theme_bw() +
      ylim(-0.15, NA)
  }else{
    ggdendrogram(dendrogramme(), label = FALSE) +
      geom_text(data = info_dendrogramme(),
                aes_string(label = colnames(info_dendrogramme()[3]), x = colnames(info_dendrogramme()[1]), y = 0, colour = input$selectVariable),
                angle = 75,
                vjust = 0,
                nudge_y = -0.1,
                size = 3) +
      theme_bw() +
      scale_colour_gradientn(colours = terrain.colors(10)) +
      ylim(-0.15, NA)
  }
  
})

# PCoA
output$PCoA <- renderPlotly({
  
  if (is.numeric(sampleAnnot_2()[,input$selectVariable3] )){
    ggplot(data = PCoA(),
           aes_string(x = colnames(PCoA())[1], y = colnames(PCoA())[2], col = input$selectVariable3)) +
      geom_point(size = 4) +
      theme_bw() +
      scale_colour_gradientn(colours = terrain.colors(10)) +
      labs(x = sprintf("Axis 1 [%s%% Variance]", 100 * round(dudi_PCoA()$eig[1] / sum(dudi_PCoA()$eig), 2)),
           y = sprintf("Axis 2 [%s%% Variance]", 100 * round(dudi_PCoA()$eig[2] / sum(dudi_PCoA()$eig), 2)))
  }else{
    
    ggplot(data = PCoA(),
           aes_string(x = colnames(PCoA())[1], y = colnames(PCoA())[2], col = input$selectVariable3)) +
      geom_point(size = 4) +
      theme_bw() +
      labs(x = sprintf("Axis 1 [%s%% Variance]", 100 * round(dudi_PCoA()$eig[1] / sum(dudi_PCoA()$eig), 2)),
           y = sprintf("Axis 2 [%s%% Variance]", 100 * round(dudi_PCoA()$eig[2] / sum(dudi_PCoA()$eig), 2)))
  }
  
})


# General Relative Abundance Dataset 2
output$Relative_Abundance_sec <- renderPlot({
  ggplot(exprDat_Relative_abundance_sec(), aes_string(x = "variable", y = "value", fill = input$selectTaxo_sec)) +
    geom_bar(stat = "identity") +
    labs(x = "Samples", 
         y = "Relative Abundance", 
         title = paste("Relative abundance at the ", input$selectTaxo_sec, " Taxonomic level", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
})

# Dendrogramme Dataset 2
output$Dendrogramme_sec <- renderPlot({
  if (!is.numeric(sampleAnnot_2()[,input$selectVariable_sec])){
    ggdendrogram(dendrogramme_sec(), label = FALSE) +
      theme_bw() +
      geom_text(data = info_dendrogramme_sec(),
                aes_string(label = colnames(info_dendrogramme_sec()[3]), x = colnames(info_dendrogramme_sec()[1]), y = 0, colour = input$selectVariable_sec),
                angle = 75,
                vjust = 0,
                nudge_y = -0.1,
                size = 3) +
      ylim(-0.15, NA)
  }else{
    ggdendrogram(dendrogramme_sec(), label = FALSE) +
      geom_text(data = info_dendrogramme_sec(),
                aes_string(label = colnames(info_dendrogramme_sec()[3]), x = colnames(info_dendrogramme_sec()[1]), y = 0, colour = input$selectVariable_sec),
                angle = 75,
                vjust = 0,
                nudge_y = -0.1,
                size = 3) +
      theme_bw() +
      scale_colour_gradientn(colours = terrain.colors(10)) +
      ylim(-0.15, NA)
  }
  
})

# PCoA Dataset 2
output$PCoA_sec <- renderPlotly({
  if (!is.numeric(sampleAnnot_2()[,input$selectVariable3_sec])){
    ggplot(data = PCoA_sec(),
           aes_string(x = colnames(PCoA_sec())[1], y = colnames(PCoA_sec())[2], col = input$selectVariable3_sec)) +
      geom_point(size = 4) +
      theme_bw() +
      labs(x = sprintf("Axis 1 [%s%% Variance]", 100 * round(dudi_PCoA_sec()$eig[1] / sum(dudi_PCoA_sec()$eig), 2)),
           y = sprintf("Axis 2 [%s%% Variance]", 100 * round(dudi_PCoA_sec()$eig[2] / sum(dudi_PCoA_sec()$eig), 2)))
  }else{
    ggplot(data = PCoA_sec(),
           aes_string(x = colnames(PCoA_sec())[1], y = colnames(PCoA_sec())[2], col = input$selectVariable3_sec)) +
      geom_point(size = 4) +
      theme_bw() +
      scale_colour_gradientn(colours = terrain.colors(10)) +
      labs(x = sprintf("Axis 1 [%s%% Variance]", 100 * round(dudi_PCoA_sec()$eig[1] / sum(dudi_PCoA_sec()$eig), 2)),
           y = sprintf("Axis 2 [%s%% Variance]", 100 * round(dudi_PCoA_sec()$eig[2] / sum(dudi_PCoA_sec()$eig), 2)))
  }
})


#### DOWNLOADS ####

output$Download_Page2_dataset1_rel_ab <- downloadHandler(
  filename = function() {
    paste("Rel_ab_Dataset1.zip")
  },
  content = function(filename) {
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    
    if (input$pdf_or_svg_page2_dataset1_rel_ab == "pdf"){
      if (input$TaxonFile || input$LoadExample == "Yes" || input$LoadExample2 == "Yes"){
        for (i in 1:ncol(taxTable_present())){
          pdf(paste("rel_abundance_", colnames(taxTable_present())[i],".pdf", sep= ""), width = 10, height = 10)
          print(ggplot(exprDat_Relative_abundance(), aes_string(x = "variable", y = "value", fill = colnames(taxTable_present())[i])) +
                  geom_bar(stat = "identity") +
                  labs(x = "Samples", y = "Relative Abundance",
                       title = paste("Relative abundance at the ", colnames(taxTable_present())[i], " Taxonomic level", sep = "")) +
                  theme(axis.text.x = element_text(angle = 90, hjust = 1)))
          dev.off()
          fs <- c(fs, paste("rel_abundance_", colnames(taxTable_present())[i],".pdf", sep= ""))
        }
      }
    }else{
      if (input$TaxonFile || input$LoadExample == "Yes" || input$LoadExample2 == "Yes"){
        for (i in 1:ncol(taxTable_present())){
          svg(paste("rel_abundance_", colnames(taxTable_present())[i],".svg", sep= ""), width = 10, height = 10)
          print(ggplot(exprDat_Relative_abundance(), aes_string(x = "variable", y = "value", fill = colnames(taxTable_present())[i])) +
                  geom_bar(stat = "identity") +
                  labs(x = "Samples", y = "Relative Abundance",
                       title = paste("Relative abundance at the ", colnames(taxTable_present())[i], " Taxonomic level", sep = "")) +
                  theme(axis.text.x = element_text(angle = 90, hjust = 1)))
          dev.off()
          fs <- c(fs, paste("rel_abundance_", colnames(taxTable_present())[i],".svg", sep= ""))
        }
      }
    }
    zip(zipfile=filename, files=fs)
  },
  contentType = "application/zip"
)

output$Download_Page2_dataset1_dendrogramme <- downloadHandler(
  filename = function() {
    paste("dendrogramme_Dataset1.", input$pdf_or_svg_page2_dataset1_dendrogramme, sep = "")
  },
  content = function(file) {
    if (input$pdf_or_svg_page2_dataset1_dendrogramme == "pdf"){
      pdf(file, width = 10, height = 10)
      print(ggdendrogram(dendrogramme(), label = FALSE) +
              geom_text(data = info_dendrogramme(),
                        aes_string(label = colnames(info_dendrogramme()[3]),
                                   x = colnames(info_dendrogramme()[1]),
                                   y = 0, colour = input$selectVariable),
                        angle = 75,
                        vjust = 0,
                        nudge_y = -0.1,
                        size = 3) +
              ylim(-0.15, NA))
      dev.off()
    }else{
      svg(file, width = 10, height = 10)
      print(ggdendrogram(dendrogramme(), label = FALSE) +
              geom_text(data = info_dendrogramme(),
                        aes_string(label = colnames(info_dendrogramme()[3]),
                                   x = colnames(info_dendrogramme()[1]),
                                   y = 0, colour = input$selectVariable),
                        angle = 75,
                        vjust = 0,
                        nudge_y = -0.1,
                        size = 3) +
              ylim(-0.15, NA))
      dev.off()
    }
  },
  contentType = paste("application/", input$pdf_or_svg_page2_dataset1_dendrogramme, sep= "")
)

output$Download_Page2_dataset1_ordination <- downloadHandler(
  filename = function() {
    paste("ordination_Dataset1.", input$pdf_or_svg_page2_dataset1_ordination, sep = "")
  },
  content = function(file) {
    if (input$pdf_or_svg_page2_dataset1_ordination == "pdf"){
      pdf(file, width = 10, height = 10)
      print(ggplot(data = PCoA(), aes_string(x = colnames(PCoA())[1], y = colnames(PCoA())[2], col = input$selectVariable3)) +
              geom_point(size = 2) +
              labs(x = sprintf("Axis 1 [%s%% Variance]", 100 * round(dudi_PCoA()$eig[1] / sum(dudi_PCoA()$eig), 2)),
                   y = sprintf("Axis 2 [%s%% Variance]", 100 * round(dudi_PCoA()$eig[2] / sum(dudi_PCoA()$eig), 2))))
      dev.off()
    }else{
      svg(file, width = 10, height = 10)
      print(ggplot(data = PCoA(), aes_string(x = colnames(PCoA())[1], y = colnames(PCoA())[2], col = input$selectVariable3)) +
              geom_point(size = 2) +
              labs(x = sprintf("Axis 1 [%s%% Variance]", 100 * round(dudi_PCoA()$eig[1] / sum(dudi_PCoA()$eig), 2)),
                   y = sprintf("Axis 2 [%s%% Variance]", 100 * round(dudi_PCoA()$eig[2] / sum(dudi_PCoA()$eig), 2))))
      dev.off()
    }
  },
  contentType = paste("application/", input$pdf_or_svg_page2_dataset1_ordination, sep= "")
)

output$Download_Page2_dataset2_rel_ab <- downloadHandler(
  filename = function() {
    paste("Rel_ab_Dataset2.zip")
  },
  content = function(filename) {
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    
    if (input$pdf_or_svg_page2_dataset2_rel_ab == "pdf"){
      if (input$TaxonFile || input$LoadExample == "Yes" || input$LoadExample2 == "Yes"){
        for (i in 1:ncol(taxTable_present())){
          pdf(paste("rel_abundance_", colnames(taxTable_present())[i],".pdf", sep= ""), width = 10, height = 10)
          print(ggplot(exprDat_Relative_abundance_sec(), aes_string(x = "variable", y = "value", fill = colnames(taxTable_present())[i])) +
                  geom_bar(stat = "identity") +
                  labs(x = "Samples", y = "Relative Abundance",
                       title = paste("Relative abundance at the ", colnames(taxTable_present())[i], " Taxonomic level", sep = "")) +
                  theme(axis.text.x = element_text(angle = 90, hjust = 1)))
          dev.off()
          fs <- c(fs, paste("rel_abundance_", colnames(taxTable_present())[i],".pdf", sep= ""))
        }
      }
    }else{
      if (input$TaxonFile || input$LoadExample == "Yes" || input$LoadExample2 == "Yes"){
        for (i in 1:ncol(taxTable_present())){
          svg(paste("rel_abundance_", colnames(taxTable_present())[i],".svg", sep= ""), width = 10, height = 10)
          print(ggplot(exprDat_Relative_abundance_sec(), aes_string(x = "variable", y = "value", fill = colnames(taxTable_present())[i])) +
                  geom_bar(stat = "identity") +
                  labs(x = "Samples", y = "Relative Abundance",
                       title = paste("Relative abundance at the ", colnames(taxTable_present())[i], " Taxonomic level", sep = "")) +
                  theme(axis.text.x = element_text(angle = 90, hjust = 1)))
          dev.off()
          fs <- c(fs, paste("rel_abundance_", colnames(taxTable_present())[i],".svg", sep= ""))
        }
      }
    }
    zip(zipfile=filename, files=fs)
  },
  contentType = "application/zip"
)

output$Download_Page2_dataset2_dendrogramme <- downloadHandler(
  filename = function() {
    paste("dendrogramme_Dataset2.", input$pdf_or_svg_page2_dataset2_dendrogramme, sep = "")
  },
  content = function(file) {
    if (input$pdf_or_svg_page2_dataset2_dendrogramme == "pdf"){
      pdf(file, width = 10, height = 10)
      print(ggdendrogram(dendrogramme_sec(), label = FALSE) +
              geom_text(data = info_dendrogramme_sec(),
                        aes_string(label = colnames(info_dendrogramme_sec()[3]), x = colnames(info_dendrogramme_sec()[1]), y = 0, colour = input$selectVariable_sec),
                        angle = 75,
                        vjust = 0,
                        nudge_y = -0.1,
                        size = 3) +
              ylim(-0.15, NA))
      dev.off()
    }else{
      svg(file, width = 10, height = 10)
      print(ggdendrogram(dendrogramme_sec(), label = FALSE) +
              geom_text(data = info_dendrogramme_sec(),
                        aes_string(label = colnames(info_dendrogramme_sec()[3]), x = colnames(info_dendrogramme_sec()[1]), y = 0, colour = input$selectVariable_sec),
                        angle = 75,
                        vjust = 0,
                        nudge_y = -0.1,
                        size = 3) +
              ylim(-0.15, NA))
      dev.off()
    }
  },
  contentType = paste("application/", input$pdf_or_svg_page2_dataset2_dendrogramme, sep= "")
)

output$Download_Page2_dataset2_ordination <- downloadHandler(
  filename = function() {
    paste("ordination_Dataset2.", input$pdf_or_svg_page2_dataset2_ordination, sep = "")
  },
  content = function(file) {
    if (input$pdf_or_svg_page2_dataset2_ordination == "pdf"){
      pdf(file, width = 10, height = 10)
      print(ggplot(data = PCoA_sec(),
                   aes_string(x = colnames(PCoA_sec())[1], y = colnames(PCoA_sec())[2], col = input$selectVariable3_sec)) +
              geom_point(size = 2) +
              labs(x = sprintf("Axis 1 [%s%% Variance]", 100 * round(dudi_PCoA_sec()$eig[1] / sum(dudi_PCoA_sec()$eig), 2)),
                   y = sprintf("Axis 2 [%s%% Variance]", 100 * round(dudi_PCoA_sec()$eig[2] / sum(dudi_PCoA_sec()$eig), 2))))
      dev.off()
    }else{
      svg(file, width = 10, height = 10)
      print(ggplot(data = PCoA_sec(),
                   aes_string(x = colnames(PCoA_sec())[1], y = colnames(PCoA_sec())[2], col = input$selectVariable3_sec)) +
              geom_point(size = 2) +
              labs(x = sprintf("Axis 1 [%s%% Variance]", 100 * round(dudi_PCoA_sec()$eig[1] / sum(dudi_PCoA_sec()$eig), 2)),
                   y = sprintf("Axis 2 [%s%% Variance]", 100 * round(dudi_PCoA_sec()$eig[2] / sum(dudi_PCoA_sec()$eig), 2))))
      dev.off()
    }
  },
  contentType = paste("application/", input$pdf_or_svg_page2_dataset2_ordination, sep= "")
)
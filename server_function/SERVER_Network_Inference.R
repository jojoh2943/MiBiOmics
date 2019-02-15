##############################################
#### SERVER FUNCTIONS : NETWORK INFERENCE ####
##############################################
# This script contains all the function related to the network inference with WGCNA in the third tabulation of 
# MiBiOmics called network inference. The function are used to change the parameter of the network (the soft 
# power, the minimum module size, the type of correlation.)

#### REACTIVE OBJECTS ####

# Choose a set of soft-thresholding powers
sft <- reactive({
  if (ncol(exprDat_WGCNA()) > 2500){
    validate(
      need(input$goPower != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 1 button and wait. The operation may take a while.')
    )
    if (input$goPower == 0){
      return(0)
    }
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    # Call the network topology analysis function
    sft <- isolate({
      withProgress({
        pickSoftThreshold(exprDat_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed)
      }, message = 'Compute Soft Power', value = 0, detail = 'This may take a while')
    })
  }else{
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    sft <- pickSoftThreshold(exprDat_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed)
  }
  sft
})


selectedDissTOM <- reactive({
  if (ncol(exprDat_WGCNA()) > 2500){
    validate(
      need(input$goNet != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 1 button and wait. The operation may take a while.')
    )
    if (input$goNet == 0){
      return(0)
    }
    
    dissTOM <- isolate({ 
      withProgress({
        softPower <- as.numeric(input$selectPower)
        if (input$corNetwork == 'bicor'){
          adjacency <- adjacency(exprDat_WGCNA(), power = softPower, type = input$signed, corFnc = 'bicor', corOptions = list(use = 'p', maxPOutliers =0.1))
        }else{
          if (input$corNetwork == 'sparCC'){
            adjacency <- X_sparcc()
          }else{
            adjacency <- adjacency(exprDat_WGCNA(), power = softPower, type = input$signed, corOptions= paste("method = '", input$corNetwork,"'", sep = ""))
            
          }
        }
        TOM = TOMsimilarity(adjacency, TOMType = input$signed)
        1-TOM
      }, message = 'Creating Adjacency Matrix', value = 0, detail = 'This may take a while')
    })
    
  }else{
    softPower <- as.numeric(input$selectPower)
    
    if (input$corNetwork == 'bicor'){
      adjacency <- adjacency(exprDat_WGCNA(), power = softPower, type = input$signed, corFnc = 'bicor', corOptions = list(maxPOutliers =0.1))
    }else{
      if (input$corNetwork == 'sparCC'){
        adjacency <- X_sparcc()
      }else{
        adjacency <- adjacency(exprDat_WGCNA(), power = softPower, type = input$signed, corOptions= paste("method = '", input$corNetwork,"'", sep = ""))
        
      }
    }
    TOM = TOMsimilarity(adjacency)
    dissTOM = 1-TOM
  }
  
  dissTOM
})

multiScaling <- reactive({
  
  if (ncol(exprDat_WGCNA()) > 2500){
    validate(
      need(input$goNet != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
    if (input$goNet == 0){
      return(0)
    }
    
    cmd1 <- isolate({ 
      withProgress({
        cmdscale(as.dist(selectedDissTOM()),3)
      }, message = 'Multiscaling ...', value = 0, detail = 'This may take a while')
    })
  }else{
    # We also propose to use classical multi-dimensional scaling plots for visualizing the network. 
    #Here we chose 3 scaling dimensions
    cmd1=cmdscale(as.dist(selectedDissTOM()),3)
    cmd1 
  }
  
})

selectedTree <- reactive({
  if (ncol(exprDat_WGCNA()) > 2500){
    if (input$goNet == 0){
      return(0)
    }
    geneTree <- isolate({
      withProgress(message = 'Clustering...', value = 0, {
        hclust(as.dist(selectedDissTOM()), method = "average")
      })
    })
  }else{
    geneTree = hclust(as.dist(selectedDissTOM()), method = "average")
  }
  geneTree
})

selectedDynamicColor <- reactive({
  minModuleSize = input$selectModuleSize
  if (ncol(exprDat_WGCNA()) > 2500){
    validate(
      need(input$goNet != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
    if (input$goNet == 0){
      return(0)
    }
    dynamicColors <- isolate({
      withProgress(message = 'Assign Color Module...', value = 0, {
        dynamicMods <- cutreeDynamic(dendro = selectedTree(), distM = selectedDissTOM(),
                                     deepSplit = 2, pamRespectsDendro = FALSE,
                                     minClusterSize = minModuleSize)
        labels2colors(dynamicMods)
      })
    })
  }else{
    dynamicMods <- cutreeDynamic(dendro = selectedTree(), distM = selectedDissTOM(),
                                 deepSplit = 2, pamRespectsDendro = FALSE,
                                 minClusterSize = minModuleSize)
    dynamicColors = labels2colors(dynamicMods)
  }
  dynamicColors
})


selectedMEs <- reactive ({
  if (ncol(exprDat_WGCNA()) > 2500){
    validate(
      need(input$goNet != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
    if (input$goNet == 0){
      return(0)
    }
    isolate({
      withProgress(message = 'Module Eigenvalues...', value = 0, {
        MEList <- moduleEigengenes(exprDat_WGCNA(), colors = selectedDynamicColor(), excludeGrey = TRUE)
        MEs = MEList$eigengenes
      })
    })
  }else{
    MEList = moduleEigengenes(exprDat_WGCNA(), colors = selectedDynamicColor(), excludeGrey = TRUE)
    MEs = MEList$eigengenes
  }
  MEs
})

selectedMETree<- reactive({
  
  validate(
    need((ncol(selectedMEs())> 1), "Not enough modules")
  )
  if (ncol(exprDat_WGCNA()) > 2500){
    validate(
      need(input$goNet != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
    if (input$goNet == 0){
      return(0)
    }
    METree <- isolate({
      withProgress(message = 'Clustering Modules...', value = 0, {
        MEDiss = 1-cor(selectedMEs())
        hclust(as.dist(MEDiss), method = "average") 
      })
    })
  }else{
    MEDiss = 1-cor(selectedMEs())
    METree = hclust(as.dist(MEDiss), method = "average") 
  }
  METree
})




# Choose a set of soft-thresholding powers
sft2 <- reactive({
  if (ncol(exprDatSec_WGCNA()) > 2500){
    validate(
      need(input$goPower2 != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 1 button and wait. The operation may take a while.')
    )
    if (input$goPower2 == 0){
      return(0)
    }
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    # Call the network topology analysis function
    sft <- isolate({
      withProgress({
        pickSoftThreshold(exprDatSec_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed2)
      }, message = 'Compute Soft Power', value = 0, detail = 'This may take a while')
    })
  }else{
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    sft <- pickSoftThreshold(exprDatSec_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed2)
  }
  sft
})


selectedDissTOM2 <- reactive({
  if (ncol(exprDatSec_WGCNA()) > 2500){
    if (input$goNet2 == 0){
      return(0)
    }
    
    dissTOM <- isolate({ 
      withProgress({
        softPower <- as.numeric(input$selectPower2)
        if (input$corNetwork2 == 'bicor'){
          adjacency <- adjacency(exprDatSec_WGCNA(), power = softPower, type = input$signed2, corFnc = 'bicor', corOptions = list(use = 'p', maxPOutliers =0.1))
        }else{
          if (input$corNetwork2 == 'sparCC'){
            adjacency <- X_sparcc()
          }else{
            adjacency <- adjacency(exprDatSec_WGCNA(), power = softPower, type = input$signed2, corOptions= paste("method = '", input$corNetwork2,"'", sep = ""))
            
          }
        }
        TOM = TOMsimilarity(adjacency, TOMType = input$signed2)
        1-TOM
      }, message = 'Creating Adjacency Matrix', value = 0, detail = 'This may take a while')
    })
    
  }else{
    softPower <- as.numeric(input$selectPower2)
    
    if (input$corNetwork2 == 'bicor'){
      adjacency <- adjacency(exprDatSec_WGCNA(), power = softPower, type = input$signed2, corFnc = 'bicor', corOptions = list(maxPOutliers =0.1))
    }else{
      if (input$corNetwork2 == 'sparCC'){
        adjacency <- X_sparcc()
      }else{
        adjacency <- adjacency(exprDatSec_WGCNA(), power = softPower, type = input$signed2, corOptions= paste("method = '", input$corNetwork2,"'", sep = ""))
        
      }
    }
    TOM = TOMsimilarity(adjacency)
    dissTOM = 1-TOM
  }
  
  dissTOM
})

multiScaling2 <- reactive({
  # We also propose to use classical multi-dimensional scaling plots for visualizing the network. 
  #Here we chose 3 scaling dimensions
  cmd1=cmdscale(as.dist(selectedDissTOM2()),3)
  cmd1
})

selectedTree2 <- reactive({
  if (ncol(exprDatSec_WGCNA()) > 2500){
    if (input$goNet2 == 0){
      return(0)
    }
    geneTree <- isolate({
      withProgress(message = 'Clustering...', value = 0, {
        hclust(as.dist(selectedDissTOM2()), method = "average")
      })
    })
  }else{
    geneTree = hclust(as.dist(selectedDissTOM2()), method = "average")
  }
  geneTree
})

selectedDynamicColor2 <- reactive({
  minModuleSize = input$selectModuleSize2
  if (ncol(exprDatSec_WGCNA()) > 2500){
    if (input$goNet2 == 0){
      return(0)
    }
    dynamicColors <- isolate({
      withProgress(message = 'Assign Color Module...', value = 0, {
        dynamicMods <- cutreeDynamic(dendro = selectedTree2(), distM = selectedDissTOM2(),
                                     deepSplit = 2, pamRespectsDendro = FALSE,
                                     minClusterSize = minModuleSize)
        labels2colors(dynamicMods)
      })
    })
  }else{
    dynamicMods <- cutreeDynamic(dendro = selectedTree2(), distM = selectedDissTOM2(),
                                 deepSplit = 2, pamRespectsDendro = FALSE,
                                 minClusterSize = minModuleSize)
    dynamicColors = labels2colors(dynamicMods)
  }
  dynamicColors
})


selectedMEs2 <- reactive ({
  if (ncol(exprDatSec_WGCNA()) > 2500){
    if (input$goNet2 == 0){
      return(0)
    }
    isolate({
      withProgress(message = 'Module Eigenvalues...', value = 0, {
        MEList <- moduleEigengenes(exprDatSec_WGCNA(), colors = selectedDynamicColor2(), excludeGrey = TRUE)
        MEList$eigengenes
      })
    })
  }else{
    MEList = moduleEigengenes(exprDatSec_WGCNA(), colors = selectedDynamicColor2(), excludeGrey = TRUE)
    MEs = MEList$eigengenes
  }
  MEs
})

selectedMETree2 <- reactive({
  validate(
    need((ncol(selectedMEs2())> 1), "Not enough modules")
  )
  if (ncol(exprDatSec_WGCNA()) > 2500){
    if (input$goNet2 == 0){
      return(0)
    }
    METree <- isolate({
      withProgress(message = 'Clustering Modules...', value = 0, {
        MEDiss = 1-cor(selectedMEs2())
        hclust(as.dist(MEDiss), method = "average") 
      })
    })
  }else{
    MEDiss = 1-cor(selectedMEs2())
    METree = hclust(as.dist(MEDiss), method = "average") 
  }
  METree
})

#### OBSERVE EVENTS ####

observeEvent(input$showPower, {
  showModal(modalDialog(
    title = "How to choose the right Soft Power ?",
    HTML(" You need to choose the lowest power for which the scale-free topology fit index curve flattens out upon reaching a high value (generally 0.90). If the scale-free topology fit index fails to reach values above 0.8 for reasonable powers (less than 15), chances are that the data exhibit a strong driver that makes a subset of the samples globally different from the rest. The difference causes high correlation among large groups of genes which invalidates the assumption of the scale-free topology approximation."),
    easyClose = TRUE,
    size = "l"
  ))
})

observeEvent(input$showPower2, {
  showModal(modalDialog(
    title = "How to choose the right Soft Power ?",
    HTML(" You need to choose the lowest power for which the scale-free topology fit index curve flattens out upon reaching a high value (generally 0.90). If the scale-free topology fit index fails to reach values above 0.8 for reasonable powers (less than 15), chances are that the data exhibit a strong driver that makes a subset of the samples globally different from the rest. The difference causes high correlation among large groups of genes which invalidates the assumption of the scale-free topology approximation."),
    easyClose = TRUE,
    size = "l"
  ))
})


observeEvent(input$showPower_P6, {
  showModal(modalDialog(
    title = "How to choose the right Soft Power ?",
    HTML(" You need to choose the lowest power for which the scale-free topology fit index curve flattens out upon reaching a high value (generally 0.90). If the scale-free topology fit index fails to reach values above 0.8 for reasonable powers (less than 15), chances are that the data exhibit a strong driver that makes a subset of the samples globally different from the rest. The difference causes high correlation among large groups of genes which invalidates the assumption of the scale-free topology approximation."),
    easyClose = TRUE,
    size = "l"
  ))
})

observeEvent(input$showPower_D2_P6, {
  showModal(modalDialog(
    title = "How to choose the right Soft Power ?",
    HTML(" You need to choose the lowest power for which the scale-free topology fit index curve flattens out upon reaching a high value (generally 0.90). If the scale-free topology fit index fails to reach values above 0.8 for reasonable powers (less than 15), chances are that the data exhibit a strong driver that makes a subset of the samples globally different from the rest. The difference causes high correlation among large groups of genes which invalidates the assumption of the scale-free topology approximation."),
    easyClose = TRUE,
    size = "l"
  ))
})

observeEvent(input$showSigned, {
  showModal(modalDialog(
    title = "What is the difference between a signed and an unsigned network ?",
    HTML("The base of a netwotk is its adjacency matrix. In WGCNA, the adjacency matrix is computed using a co-expression similarity. 
         This quantitative value is strongly linked to the correlation between two genes profiles. In a signed network, the co-expression
         value is defined to according to the sign of the relationship between the two genes. In an unsigned network, it is the opposite. <br> <br>"),
    HTML("For unsigned Network : adjacency = |cor|^power <br> <br> "),
    HTML("For signed Netwotk : adjacency = (0.5*(1+cor))^power <br> <br> "),
    HTML("Generally, signed Network are advised for biological analysis because they tend to be closer to the system reality."),
    easyClose = TRUE,
    size = "l"
    ))
})

observeEvent(input$showSigned2, {
  showModal(modalDialog(
    title = "What is the difference between a signed and an unsigned network ?",
    HTML("The base of a netwotk is its adjacency matrix. In WGCNA, the adjacency matrix is computed using a co-expression similarity. 
         This quantitative value is strongly linked to the correlation between two genes profiles. In a signed network, the co-expression
         value is defined to according to the sign of the relationship between the two genes. In an unsigned network, it is the opposite. <br> <br>"),
    HTML("For unsigned Network : adjacency = |cor|^power <br> <br> "),
    HTML("For signed Netwotk : adjacency = (0.5*(1+cor))^power <br> <br> "),
    HTML("Generally, signed Network are advised for biological analysis because they tend to be closer to the system reality."),
    easyClose = TRUE,
    size = "l"
    ))
})

observeEvent(input$showSigned_P6, {
  showModal(modalDialog(
    title = "What is the difference between a signed and an unsigned network ?",
    HTML("The base of a netwotk is its adjacency matrix. In WGCNA, the adjacency matrix is computed using a co-expression similarity. 
         This quantitative value is strongly linked to the correlation between two genes profiles. In a signed network, the co-expression
         value is defined to according to the sign of the relationship between the two genes. In an unsigned network, it is the opposite. <br> <br>"),
    HTML("For unsigned Network : adjacency = |cor|^power <br> <br> "),
    HTML("For signed Netwotk : adjacency = (0.5*(1+cor))^power <br> <br> "),
    HTML("Generally, signed Network are advised for biological analysis because they tend to be closer to the system reality."),
    easyClose = TRUE,
    size = "l"
    ))
})

observeEvent(input$showSigned_D2_P6, {
  showModal(modalDialog(
    title = "What is the difference between a signed and an unsigned network ?",
    HTML("The base of a netwotk is its adjacency matrix. In WGCNA, the adjacency matrix is computed using a co-expression similarity. 
           This quantitative value is strongly linked to the correlation between two genes profiles. In a signed network, the co-expression
           value is defined to according to the sign of the relationship between the two genes. In an unsigned network, it is the opposite. <br> <br>"),
    HTML("For unsigned Network : adjacency = |cor|^power <br> <br> "),
    HTML("For signed Netwotk : adjacency = (0.5*(1+cor))^power <br> <br> "),
    HTML("Generally, signed Network are advised for biological analysis because they tend to be closer to the system reality."),
    easyClose = TRUE,
    size = "l"
  ))
})

observeEvent(input$Download_Data_Exploration, {
  updateTabsetPanel(session, "Page2_dataset1", selected = "Page2_dataset1_ordination")
  updateTabsetPanel(session, "Page2_dataset1", selected = "Page2_dataset1_dendrogramme")
  if (input$TypeAnalysis == "multivariate"){
    updateTabsetPanel(session, "Page2_datasets", selected = "Page2_dataset2")
    updateTabsetPanel(session, "Page2_dataset2", selected = "Page2_dataset2_ordination")
    updateTabsetPanel(session, "Page2_dataset2", selected = "Page2_dataset2_dendrogramme")
  }
  
})

#### INTERFACE VARIABLES ####


output$Power <- renderUI({
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  selectInput("selectPower", 
              label = "Choose the soft power: ", 
              choices = as.character(powers), 
              selected = as.character(sft()[["powerEstimate"]]))
})

output$Power2 <- renderUI({
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  selectInput("selectPower2", 
              label = "Choose the soft power: ", 
              choices = as.character(powers), 
              selected = as.character(sft2()[["powerEstimate"]]))
})

output$Power_P6 <- renderUI({
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  selectInput("Power_P6", 
              label = "Choose the soft power (first dataset): ", 
              choices = as.character(powers), 
              selected = as.character(sft()[["powerEstimate"]]))
})

output$Power_D2_P6 <- renderUI({
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  selectInput("Power_D2_P6", 
              label = "Choose the soft power (second dataset): ", 
              choices = as.character(powers), 
              selected = as.character(sft2()[["powerEstimate"]]))
})

output$ModuleSize <- renderUI({
  numericInput("selectModuleSize", 
               label = "Choose the minimum module size :",
               value = 6, 
               min = 2,
               max = max(30, round(ncol(exprDat_WGCNA())/3))
  )
})

output$ModuleSize2 <- renderUI({
  numericInput("selectModuleSize2", 
               label = "Choose the minimum module size :",
               value = 6, 
               min = 2,
               max = max(30, round(ncol(exprDatSec_WGCNA())/3))
  )
})

output$ModuleSize_P6 <- renderUI({
  numericInput("ModuleSize_P6", 
               label = "Choose the minimum module size (first dataset):",
               value = 6, 
               min = 2,
               max = max(30, round(ncol(exprDat_WGCNA())/3))
  )
})

output$ModuleSize_D2_P6 <- renderUI({
  numericInput("ModuleSize_D2_P6", 
               label = "Choose the minimum module size (Second dataset):",
               value = 6, 
               min = 2,
               max = max(30, round(ncol(exprDatSec_WGCNA())/3))
  )
})

# A button to trigger the net creation (only when large datasets are uploaded)
output$buttonPower <- renderUI({
  if (ncol(exprDat_WGCNA()) > 2500){
    actionButton("goPower", "STEP 1: Compute Soft Power Calculation !")
  }
})

# A button to trigger the net creation (only when large datasets are uploaded)
output$buttonPower2 <- renderUI({
  if (ncol(exprDat_WGCNA()) > 2500){
    actionButton("goPower2", "STEP 1: Compute Soft Power Calculation !")
  }
})

output$buttonNet <- renderUI({
  if (ncol(exprDat_WGCNA()) > 2500){
    actionButton("goNet", "STEP 2: Begin Network Creation !")
  }
})

output$buttonNet2 <- renderUI({
  if (ncol(exprDatSec_WGCNA()) > 2500){
    actionButton("goNet2", "STEP 3: Begin Network Creation !")
    
  }
})


#### PLOT OUTPUTS ####

#Scale Free Topology
output$ScaleFree <- renderPlot({
  ggplot(data =sft()$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq )) + 
    geom_text(aes(label = Power, color = "tomato3"), show.legend =  FALSE) +
    geom_hline(yintercept = 0.9) +
    theme_bw() +
    labs(y = "Scale Free Topology", x = "Soft Threshold (Power)", title = "Scale Independence") +
    theme(legend.position = "none")
})

output$ScaleFree_P6 <- renderPlot({
  ggplot(data =sft()$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq )) + 
    geom_text(aes(label = Power, color = "tomato3"), show.legend =  FALSE) +
    geom_hline(yintercept = 0.9) +
    theme_bw() +
    labs(y = "Scale Free Topology", x = "Soft Threshold (Power)", title = "Scale Independence") +
    theme(legend.position = "none")
})

# Module Clustering
output$ModuleClustering <- renderPlot({
  ggdendrogram(selectedMETree()) + 
    labs(y = "Height", title = "Clustering of module eigengenes")
})

# Dendrogramme Genes with Modules Colors
output$Dendrogramme_Colors <- renderPlot({
  plotDendroAndColors(selectedTree(), selectedDynamicColor(), "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
})

# Module Informations
output$Module_Information <- renderDataTable({
  datatable(as.data.frame(table(selectedDynamicColor())), 
            class = 'cell-border stripe', 
            rownames = FALSE, 
            options = list(pageLength = 8, dom = 'tip'),
            colnames = c('Module Color' = 'Var1', 'Number of element' = 'Freq'))
})

output$multiScalePlot <- renderScatterplotThree({
  scatterplot3js(multiScaling()[,1],
                 multiScaling()[,2],
                 multiScaling()[,3],
                 color = selectedDynamicColor(),
                 x.ticklabs = "Scaling Axis 1",
                 y.ticklabs = "Scaling Axis 2",
                 z.ticklabs = "Scaling Axis 3",
                 size =0.4)
})

output$multiScalePlotOutput <- renderUI({
  scatterplotThreeOutput("multiScalePlot")
})


#Scale Free Topology
output$ScaleFree2 <- renderPlot({
  ggplot(data =sft2()$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq )) + 
    geom_text(aes(label = Power, color = "tomato3"), show.legend =  FALSE) +
    geom_hline(yintercept = 0.9) +
    theme_bw() +
    labs(y = "Scale Free Topology", x = "Soft Threshold (Power)", title = "Scale Independence") +
    theme(legend.position = "none")
})

output$ScaleFree_D2_P6 <- renderPlot({
  ggplot(data =sft2()$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq )) + 
    geom_text(aes(label = Power, color = "tomato3"), show.legend =  FALSE) +
    geom_hline(yintercept = 0.9) +
    theme_bw() +
    labs(y = "Scale Free Topology", x = "Soft Threshold (Power)", title = "Scale Independence") +
    theme(legend.position = "none")
})

# Module Clustering
output$ModuleClustering2 <- renderPlot({
  ggdendrogram(selectedMETree2()) + 
    labs(y = "Height", title = "Clustering of module eigengenes")
})

# Dendrogramme Genes with Modules Colors
output$Dendrogramme_Colors2 <- renderPlot({
  plotDendroAndColors(selectedTree2(), selectedDynamicColor2(), "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
})

# Module Informations
output$Module_Information2 <- renderDataTable({
  datatable(as.data.frame(table(selectedDynamicColor2())), 
            class = 'cell-border stripe', 
            rownames = FALSE, 
            options = list(pageLength = 8, dom = 'tip'),
            colnames = c('Module Color' = 'Var1', 'Number of element' = 'Freq'))
})

output$multiScalePlot2 <- renderScatterplotThree({
  scatterplot3js(multiScaling2()[,1], 
                 multiScaling2()[,2], 
                 multiScaling2()[,3], 
                 color = selectedDynamicColor2(), 
                 x.ticklabs = "Scaling Axis 1", 
                 y.ticklabs = "Scaling Axis 2", 
                 z.ticklabs = "Scaling Axis 3", 
                 size =0.4)
})

output$multiScalePlotOutput2 <- renderUI({
  scatterplotThreeOutput("multiScalePlot2")
})



#### DOWNLOADS ####

output$Download_Network_Creation <- downloadHandler(
  filename = function(){
    paste("NETWORK_CREATION_dataset1.zip")
  },
  content = function (filename){
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    if (input$pdf_or_svg_p3 == "pdf"){
      pdf("Soft_Power.pdf", width = 10, height = 10)
      print(ggplot(data =sft()$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq )) + 
              geom_text(aes(label = Power, color = "tomato3"), show.legend =  FALSE) +
              geom_hline(yintercept = 0.9) +
              labs(y = "Scale Free Topology", x = "Soft Threshold (Power)", title = "Scale Independence") +
              theme(legend.position = "none"))
      dev.off()
      
      fs <- c(fs, "Soft_Power.pdf")
      pdf("Module_dendrogramme.pdf", width = 10, height = 10)
      print(ggdendrogram(selectedMETree()) + 
              labs(y = "Height", title = "Clustering of module eigengenes"))
      dev.off()
      fs <- c(fs, "Module_dendrogramme.pdf")
      pdf("gene_dendrogramme_colors.pdf", width = 10, height = 10)
      print(plotDendroAndColors(selectedTree(), selectedDynamicColor(), "Dynamic Tree Cut",
                                dendroLabels = FALSE, hang = 0.03,
                                addGuide = TRUE, guideHang = 0.05,
                                main = "Gene dendrogram and module colors"))
      dev.off()
      fs <- c(fs, "gene_dendrogramme_colors.pdf")
    }else{
      svg("Soft_Power.svg", width = 10, height = 10)
      print(ggplot(data =sft()$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq )) + 
              geom_text(aes(label = Power, color = "tomato3"), show.legend =  FALSE) +
              geom_hline(yintercept = 0.9) +
              labs(y = "Scale Free Topology", x = "Soft Threshold (Power)", title = "Scale Independence") +
              theme(legend.position = "none"))
      dev.off()
      
      fs <- c(fs, "Soft_Power.svg")
      svg("Module_dendrogramme.svg", width = 10, height = 10)
      print(ggdendrogram(selectedMETree()) + 
              labs(y = "Height", title = "Clustering of module eigengenes"))
      dev.off()
      fs <- c(fs, "Module_dendrogramme.svg")
      svg("gene_dendrogramme_colors.svg", width = 10, height = 10)
      print(plotDendroAndColors(selectedTree(), selectedDynamicColor(), "Dynamic Tree Cut",
                                dendroLabels = FALSE, hang = 0.03,
                                addGuide = TRUE, guideHang = 0.05,
                                main = "Gene dendrogram and module colors"))
      dev.off()
      fs <- c(fs, "gene_dendrogramme_colors.svg")
    }
    
    zip(zipfile=filename, files=fs)
    
  },
  contentType = "application/zip")


output$Download_Network_Creation_dataset2 <- downloadHandler(
  filename = function(){
    paste("NETWORK_CREATION_dataset2.zip")
  },
  content = function (filename){
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    if (input$pdf_or_svg_p3_dataset2 == "pdf"){
      pdf("Soft_Power.pdf", width = 10, height = 10)
      print(ggplot(data =sft2()$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq )) + 
              geom_text(aes(label = Power, color = "tomato3"), show.legend =  FALSE) +
              geom_hline(yintercept = 0.9) +
              labs(y = "Scale Free Topology", x = "Soft Threshold (Power)", title = "Scale Independence") +
              theme(legend.position = "none"))
      dev.off()
      
      fs <- c(fs, "Soft_Power.pdf")
      pdf("Module_dendrogramme.pdf", width = 10, height = 10)
      print(ggdendrogram(selectedMETree2()) + 
              labs(y = "Height", title = "Clustering of module eigengenes"))
      dev.off()
      fs <- c(fs, "Module_dendrogramme.pdf")
      pdf("gene_dendrogramme_colors.pdf", width = 10, height = 10)
      print(plotDendroAndColors(selectedTree2(), selectedDynamicColor2(), "Dynamic Tree Cut",
                                dendroLabels = FALSE, hang = 0.03,
                                addGuide = TRUE, guideHang = 0.05,
                                main = "Gene dendrogram and module colors"))
      dev.off()
      fs <- c(fs, "gene_dendrogramme_colors.pdf")
    }else{
      svg("Soft_Power.svg", width = 10, height = 10)
      print(ggplot(data =sft2()$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq )) + 
              geom_text(aes(label = Power, color = "tomato3"), show.legend =  FALSE) +
              geom_hline(yintercept = 0.9) +
              labs(y = "Scale Free Topology", x = "Soft Threshold (Power)", title = "Scale Independence") +
              theme(legend.position = "none"))
      dev.off()
      
      fs <- c(fs, "Soft_Power.svg")
      svg("Module_dendrogramme.svg", width = 10, height = 10)
      print(ggdendrogram(selectedMETree2()) + 
              labs(y = "Height", title = "Clustering of module eigengenes"))
      dev.off()
      fs <- c(fs, "Module_dendrogramme.svg")
      svg("gene_dendrogramme_colors.svg", width = 10, height = 10)
      print(plotDendroAndColors(selectedTree2(), selectedDynamicColor2(), "Dynamic Tree Cut",
                                dendroLabels = FALSE, hang = 0.03,
                                addGuide = TRUE, guideHang = 0.05,
                                main = "Gene dendrogram and module colors"))
      dev.off()
      fs <- c(fs, "gene_dendrogramme_colors.svg")
    }
    
    zip(zipfile=filename, files=fs)
    
  },
  contentType = "application/zip")


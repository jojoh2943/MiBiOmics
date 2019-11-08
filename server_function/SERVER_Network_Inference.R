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
        if(input$corNetwork == 'bicor'){
           pickSoftThreshold(exprDat_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed, corFnc = 'bicor', corOptions = list(use = 'p', maxPOutliers =0.1))
        }else{
          
          sft <- pickSoftThreshold(exprDat_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed, corOptions = list(use = 'p', method = as.character(input$corNetwork)))
        }      }, message = 'Compute Soft Power', value = 0, detail = 'This may take a while')
    })
  }else{
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    if(input$corNetwork == 'bicor'){
      sft <- pickSoftThreshold(exprDat_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed, corFnc = 'bicor', corOptions = list(use = 'p', maxPOutliers =0.1))
    }else{
        sft <- pickSoftThreshold(exprDat_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed, corOptions = list(use = 'p', method = as.character(input$corNetwork)))

    }
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
    print(dim(adjacency))
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
    colnames(cmd1) <- c("scale1", "scale2", "scale3")
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
  rownames(MEs) <- rownames(exprDat_WGCNA())
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
        if(input$corNetwork2 == 'bicor'){
          pickSoftThreshold(exprDatSec_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed2, corFnc = 'bicor', corOptions = list(use = 'p', maxPOutliers =0.1))
        }else{
          sft <- pickSoftThreshold(exprDatSec_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed2, corOptions = list(use = 'p', method = as.character(input$corNetwork2)))
        }      
        }, message = 'Compute Soft Power', value = 0, detail = 'This may take a while')
    })
  }else{
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    if(input$corNetwork2 == 'bicor'){
      sft <- pickSoftThreshold(exprDatSec_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed2, corFnc = 'bicor', corOptions = list(use = 'p', maxPOutliers =0.1))
    }else{
      print(paste("use = 'p', method = '", input$corNetwork,"'", sep = ""))
      sft <- pickSoftThreshold(exprDatSec_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed2, corOptions = list(use = 'p', method = as.character(input$corNetwork2)))
    } 
    }
  sft
})


selectedDissTOM2 <- reactive({
  if (ncol(exprDatSec_WGCNA()) > 2500){
    validate(
      need(input$goNet2 != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
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
  
  if (ncol(exprDatSec_WGCNA()) > 2500){
    validate(
      need(input$goNet2 != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
    if (input$goNet2 == 0){
      return(0)
    }
    
    cmd1 <- isolate({ 
      withProgress({
        cmdscale(as.dist(selectedDissTOM2()),3)
      }, message = 'Multiscaling ...', value = 0, detail = 'This may take a while')
    })
  }else{
    # We also propose to use classical multi-dimensional scaling plots for visualizing the network. 
    #Here we chose 3 scaling dimensions
    cmd1=cmdscale(as.dist(selectedDissTOM2()),3)
    colnames(cmd1) <- c("scale1", "scale2", "scale3")
    cmd1
  }
})

selectedTree2 <- reactive({
  if (ncol(exprDatSec_WGCNA()) > 2500){
    validate(
      need(input$goNet2 != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
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
    validate(
      need(input$goNet2 != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
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
    validate(
      need(input$goNet2 != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
    if (input$goNet2 == 0){
      return(0)
    }
    isolate({
      withProgress(message = 'Module Eigenvalues...', value = 0, {
        MEList <- moduleEigengenes(exprDatSec_WGCNA(), colors = selectedDynamicColor2(), excludeGrey = TRUE)
        MEs <- MEList$eigengenes
      })
    })
  }else{
    MEList = moduleEigengenes(exprDatSec_WGCNA(), colors = selectedDynamicColor2(), excludeGrey = TRUE)
    MEs = MEList$eigengenes
  }
  rownames(MEs) <- rownames(exprDatSec_WGCNA())
  MEs
})

selectedMETree2 <- reactive({
  validate(
    need((ncol(selectedMEs2())> 1), "Not enough modules")
  )
  if (ncol(exprDatSec_WGCNA()) > 2500){
    validate(
      need(input$goNet2 != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
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

# Choose a set of soft-thresholding powers
sft3 <- reactive({
  if (ncol(exprDatTer_WGCNA()) > 2500){
    validate(
      need(input$goPower3 != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 1 button and wait. The operation may take a while.')
    )
    if (input$goPower3 == 0){
      return(0)
    }
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    # Call the network topology analysis function
    sft <- isolate({
      withProgress({
        if(input$corNetwork3 == 'bicor'){
          pickSoftThreshold(exprDatTer_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed3, corFnc = 'bicor', corOptions = list(use = 'p', maxPOutliers =0.1))
        }else{
          sft <- pickSoftThreshold(exprDatTer_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed3, corOptions = list(use = 'p', method = as.character(input$corNetwork3)))
        }  
        }, message = 'Compute Soft Power', value = 0, detail = 'This may take a while')
    })
  }else{
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    if(input$corNetwork3 == 'bicor'){
      sft <- pickSoftThreshold(exprDatTer_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed3, corFnc = 'bicor', corOptions = list(use = 'p', maxPOutliers =0.1))
    }else{
      sft <- pickSoftThreshold(exprDatTer_WGCNA(), powerVector = powers, verbose = 5, networkType = input$signed3, corOptions = list(use = 'p', method = as.character(input$corNetwork3)))
    } 
    }
  sft
})


selectedDissTOM3 <- reactive({
  if (ncol(exprDatTer_WGCNA()) > 2500){
    validate(
      need(input$goNet3 != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
    if (input$goNet3 == 0){
      return(0)
    }
    
    dissTOM <- isolate({ 
      withProgress({
        softPower <- as.numeric(input$selectPower3)
        if (input$corNetwork3 == 'bicor'){
          adjacency <- adjacency(exprDatTer_WGCNA(), power = softPower, type = input$signed3, corFnc = 'bicor', corOptions = list(use = 'p', maxPOutliers =0.1))
        }else{
          adjacency <- adjacency(exprDatTer_WGCNA(), power = softPower, type = input$signed3, corOptions= paste("method = '", input$corNetwork3,"'", sep = ""))
          
        }
        TOM = TOMsimilarity(adjacency, TOMType = input$signed3)
        1-TOM
      }, message = 'Creating Adjacency Matrix', value = 0, detail = 'This may take a while')
    })
    
  }else{
    softPower <- as.numeric(input$selectPower3)
    
    if (input$corNetwork3 == 'bicor'){
      adjacency <- adjacency(exprDatTer_WGCNA(), power = softPower, type = input$signed3, corFnc = 'bicor', corOptions = list(maxPOutliers =0.1))
    }else{
      adjacency <- adjacency(exprDatTer_WGCNA(), power = softPower, type = input$signed3, corOptions= paste("method = '", input$corNetwork3,"'", sep = ""))
    }
    TOM = TOMsimilarity(adjacency)
    dissTOM = 1-TOM
  }
  dissTOM
})

multiScaling3 <- reactive({
  if (ncol(exprDatTer_WGCNA()) > 2500){
    validate(
      need(input$goNet3 != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
    if (input$goNet3 == 0){
      return(0)
    }
    
    cmd1 <- isolate({ 
      withProgress({
        cmdscale(as.dist(selectedDissTOM3()),3)
      }, message = 'Multiscaling ...', value = 0, detail = 'This may take a while')
    })
  }else{
    # We also propose to use classical multi-dimensional scaling plots for visualizing the network. 
    #Here we chose 3 scaling dimensions
    cmd1=cmdscale(as.dist(selectedDissTOM3()),3)
    colnames(cmd1) <- c("scale1", "scale2", "scale3")
    cmd1
  }
})

selectedTree3 <- reactive({
  if (ncol(exprDatTer_WGCNA()) > 2500){
    validate(
      need(input$goNet3 != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
    if (input$goNet3 == 0){
      return(0)
    }
    geneTree <- isolate({
      withProgress(message = 'Clustering...', value = 0, {
        hclust(as.dist(selectedDissTOM3()), method = "average")
      })
    })
  }else{
    geneTree = hclust(as.dist(selectedDissTOM3()), method = "average")
  }
  geneTree
})

selectedDynamicColor3 <- reactive({
  minModuleSize = input$selectModuleSize3
  if (ncol(exprDatTer_WGCNA()) > 2500){
    validate(
      need(input$goNet3 != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
    if (input$goNet3 == 0){
      return(0)
    }
    dynamicColors <- isolate({
      withProgress(message = 'Assign Color Module...', value = 0, {
        dynamicMods <- cutreeDynamic(dendro = selectedTree3(), distM = selectedDissTOM3(),
                                     deepSplit = 2, pamRespectsDendro = FALSE,
                                     minClusterSize = minModuleSize)
        labels2colors(dynamicMods)
      })
    })
  }else{
    dynamicMods <- cutreeDynamic(dendro = selectedTree3(), distM = selectedDissTOM3(),
                                 deepSplit = 2, pamRespectsDendro = FALSE,
                                 minClusterSize = minModuleSize)
    dynamicColors = labels2colors(dynamicMods)
  }
  dynamicColors
})


selectedMEs3 <- reactive ({
  if (ncol(exprDatTer_WGCNA()) > 2500){
    validate(
      need(input$goNet3 != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
    if (input$goNet3 == 0){
      return(0)
    }
    isolate({
      withProgress(message = 'Module Eigenvalues...', value = 0, {
        MEList <- moduleEigengenes(exprDatTer_WGCNA(), colors = selectedDynamicColor3(), excludeGrey = TRUE)
        MEs <- MEList$eigengenes
      })
    })
  }else{
    MEList = moduleEigengenes(exprDatTer_WGCNA(), colors = selectedDynamicColor3(), excludeGrey = TRUE)
    MEs = MEList$eigengenes
  }
  rownames(MEs) <- rownames(exprDatTer_WGCNA())
  MEs
})

selectedMETree3 <- reactive({
  validate(
    need((ncol(selectedMEs3())> 1), "Not enough modules")
  )
  if (ncol(exprDatTer_WGCNA()) > 2500){
    validate(
      need(input$goNet3 != 0, 'Because your dataset is too large, we will perform each step independantly. Push the STEP 2 button and wait. The operation may take a while.')
    )
    if (input$goNet3 == 0){
      return(0)
    }
    METree <- isolate({
      withProgress(message = 'Clustering Modules...', value = 0, {
        MEDiss = 1-cor(selectedMEs3())
        hclust(as.dist(MEDiss), method = "average") 
      })
    })
  }else{
    MEDiss = 1-cor(selectedMEs3())
    METree = hclust(as.dist(MEDiss), method = "average") 
  }
  print(input$enableTabulations)
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

observeEvent(input$showPower3, {
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

observeEvent(input$showSigned3, {
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

output$Power3 <- renderUI({
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  selectInput("selectPower3", 
              label = "Choose the soft power: ", 
              choices = as.character(powers), 
              selected = as.character(sft3()[["powerEstimate"]]))
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

output$ModuleSize3 <- renderUI({
  numericInput("selectModuleSize3", 
               label = "Choose the minimum module size :",
               value = 6, 
               min = 2,
               max = max(30, round(ncol(exprDatTer_WGCNA())/3))
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
  if (ncol(exprDatSec_WGCNA()) > 2500){
    actionButton("goPower2", "STEP 1: Compute Soft Power Calculation !")
  }
})

# A button to trigger the net creation (only when large datasets are uploaded)
output$buttonPower3 <- renderUI({
  if (ncol(exprDatTer_WGCNA()) > 2500){
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
    actionButton("goNet2", "STEP 2: Begin Network Creation !")
    
  }
})

output$buttonNet3 <- renderUI({
  if (ncol(exprDatTer_WGCNA()) > 2500){
    actionButton("goNet3", "STEP 2: Begin Network Creation !")
    
  }
})


#### PLOT OUTPUTS ####

#Scale Free Topology
output$ScaleFree <- renderPlot({
  if (is.null(sft()))
    return(NULL)
  ggplot(data =sft()$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq )) + 
    geom_text(aes(label = Power, color = "tomato3"), show.legend =  FALSE) +
    geom_hline(yintercept = 0.9) +
    theme_bw() +
    labs(y = "Scale Free Topology", x = "Soft Threshold (Power)", title = "Scale Independence") +
    theme(legend.position = "none")
})


# Module Clustering
output$ModuleClustering <- renderPlot({
  if (is.null(selectedMETree()))
    return(NULL)
  ggdendrogram(selectedMETree()) + 
    labs(y = "Height", title = "Clustering of module eigengenes")
})

# Dendrogramme Genes with Modules Colors
output$Dendrogramme_Colors <- renderPlot({
  if (is.null(selectedTree()))
    return(NULL)
  plotDendroAndColors(selectedTree(), selectedDynamicColor(), "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
})

# Module Informations
output$Module_Information <- renderDataTable({
  if (is.null(selectedDynamicColor()))
    return(NULL)
  datatable(as.data.frame(table(selectedDynamicColor())), 
            class = 'cell-border stripe', 
            rownames = FALSE, 
            options = list(pageLength = 8, dom = 'tip'),
            colnames = c('Module Color' = 'Var1', 'Number of element' = 'Freq'))
})

output$multiScalePlot <- renderScatterplotThree({
  if (is.null(multiScaling()))
    return(NULL)
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
  if (is.null(multiScaling()))
    return(NULL)
  scatterplotThreeOutput("multiScalePlot")
})


#Scale Free Topology
output$ScaleFree2 <- renderPlot({
  if (is.null(sft2()))
    return(NULL)
  ggplot(data =sft2()$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq )) + 
    geom_text(aes(label = Power, color = "tomato3"), show.legend =  FALSE) +
    geom_hline(yintercept = 0.9) +
    theme_bw() +
    labs(y = "Scale Free Topology", x = "Soft Threshold (Power)", title = "Scale Independence") +
    theme(legend.position = "none")
})


# Module Clustering
output$ModuleClustering2 <- renderPlot({
  if (is.null(selectedMETree2()))
    return(NULL)
  ggdendrogram(selectedMETree2()) + 
    labs(y = "Height", title = "Clustering of module eigengenes")
})

# Dendrogramme Genes with Modules Colors
output$Dendrogramme_Colors2 <- renderPlot({
  if (is.null(selectedDynamicColor2()))
    return(NULL)
  plotDendroAndColors(selectedTree2(), selectedDynamicColor2(), "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
})

# Module Informations
output$Module_Information2 <- renderDataTable({
  if (is.null(selectedDynamicColor2()))
    return(NULL)
  datatable(as.data.frame(table(selectedDynamicColor2())), 
            class = 'cell-border stripe', 
            rownames = FALSE, 
            options = list(pageLength = 8, dom = 'tip'),
            colnames = c('Module Color' = 'Var1', 'Number of element' = 'Freq'))
})

output$multiScalePlot2 <- renderScatterplotThree({
  if (is.null(multiScaling2()))
    return(NULL)
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
  if (is.null(multiScaling2()))
    return(NULL)
  scatterplotThreeOutput("multiScalePlot2")
})

#Scale Free Topology
output$ScaleFree3 <- renderPlot({
  if (is.null(sft3()))
    return(NULL)
  ggplot(data =sft3()$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq )) + 
    geom_text(aes(label = Power, color = "tomato3"), show.legend =  FALSE) +
    geom_hline(yintercept = 0.9) +
    theme_bw() +
    labs(y = "Scale Free Topology", x = "Soft Threshold (Power)", title = "Scale Independence") +
    theme(legend.position = "none")
})


# Module Clustering
output$ModuleClustering3 <- renderPlot({
  if (is.null(selectedMETree3()))
    return(NULL)
  ggdendrogram(selectedMETree3()) + 
    labs(y = "Height", title = "Clustering of module eigengenes")
})



# Dendrogramme Genes with Modules Colors
output$Dendrogramme_Colors3 <- renderPlot({
  if (is.null(selectedDynamicColor3()))
    return(NULL)
  plotDendroAndColors(selectedTree3(), selectedDynamicColor3(), "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
})

# Module Informations
output$Module_Information3 <- renderDataTable({
  if (is.null(selectedDynamicColor3()))
    return(NULL)
  datatable(as.data.frame(table(selectedDynamicColor3())), 
            class = 'cell-border stripe', 
            rownames = FALSE, 
            options = list(pageLength = 8, dom = 'tip'),
            colnames = c('Module Color' = 'Var1', 'Number of element' = 'Freq'))
})

output$multiScalePlot3 <- renderScatterplotThree({
  if (is.null(multiScaling3()))
    return(NULL)
  scatterplot3js(multiScaling3()[,1], 
                 multiScaling3()[,2], 
                 multiScaling3()[,3], 
                 color = selectedDynamicColor3(), 
                 x.ticklabs = "Scaling Axis 1", 
                 y.ticklabs = "Scaling Axis 2", 
                 z.ticklabs = "Scaling Axis 3", 
                 size =0.4)
})

output$multiScalePlotOutput3 <- renderUI({
  if (is.null(multiScaling3()))
    return(NULL)
  scatterplotThreeOutput("multiScalePlot3")
})



#### DOWNLOADS ####


output$downloadEigen <- downloadHandler(
  filename = function(){
    "Modules_Eigenvalue.csv"
  },
  content = function (filename){
    write.csv(selectedMEs(), filename)
  }
)

output$downloadModule <- downloadHandler(
  filename = function(){
    "Modules_Size.csv"
  },
  content = function (filename){
    write.csv(as.data.frame(table(selectedDynamicColor())), filename)
  }
)

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
      pdf("Multiscaling.pdf", width = 10, height = 10)
      multiscaling <- as.data.frame(multiScaling())
      multiscaling$color <- as.vector(selectedDynamicColor())
      my_colors <- as.vector(unique(multiscaling$color))
      my_colors <- my_colors[order(my_colors)]
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale1, y = scale2, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 1",
           y = "Scaling Dimension 2")
            + scale_color_manual(values=my_colors)) 
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale2, y = scale3, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 2",
                   y = "Scaling Dimension 3")
            + scale_color_manual(values=my_colors))
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale1, y = scale3, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 1",
                   y = "Scaling Dimension 3")
            + scale_color_manual(values=my_colors))
      dev.off()

      fs <- c(fs, "Multiscaling.pdf")
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
      multiscaling <- as.data.frame(multiScaling())
      multiscaling$color <- as.vector(selectedDynamicColor())
      my_colors <- as.vector(unique(multiscaling$color))
      my_colors <- my_colors[order(my_colors)]
      svg("Multiscaling_1_2.svg", width = 10, height = 10)
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale1, y = scale2, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 1",
                   y = "Scaling Dimension 2")
            + scale_color_manual(values=my_colors)) 
      dev.off()
      
      fs <- c(fs, "Multiscaling_1_2.svg")
      
      svg("Multiscaling_2_3.svg", width = 10, height = 10)
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale2, y = scale3, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 2",
                   y = "Scaling Dimension 3")
            + scale_color_manual(values=my_colors))
      dev.off()
      
      fs <- c(fs, "Multiscaling_2_3.svg")
      
      svg("Multiscaling_1_3.svg", width = 10, height = 10)
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale1, y = scale3, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 1",
                   y = "Scaling Dimension 3")
            + scale_color_manual(values=my_colors))
      dev.off()
      
      fs <- c(fs, "Multiscaling_1_3.svg")
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

#########################"
output$downloadEigenD2 <- downloadHandler(
  filename = function(){
    "Modules_Eigenvalue.csv"
  },
  content = function (filename){
    write.csv(selectedMEs2(), filename)
  }
)

output$downloadModuleD2 <- downloadHandler(
  filename = function(){
    "Modules_Size.csv"
  },
  content = function (filename){
    write.csv(as.data.frame(table(selectedDynamicColor2())), filename)
  }
)

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
      pdf("Multiscaling.pdf", width = 10, height = 10)
      multiscaling <- as.data.frame(multiScaling2())
      multiscaling$color <- as.vector(selectedDynamicColor2())
      my_colors <- as.vector(unique(multiscaling$color))
      my_colors <- my_colors[order(my_colors)]
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale1, y = scale2, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 1",
                   y = "Scaling Dimension 2")
            + scale_color_manual(values=my_colors)) 
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale2, y = scale3, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 2",
                   y = "Scaling Dimension 3")
            + scale_color_manual(values=my_colors))
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale1, y = scale3, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 1",
                   y = "Scaling Dimension 3")
            + scale_color_manual(values=my_colors))
      dev.off()
      
      fs <- c(fs, "Multiscaling.pdf")
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
      multiscaling <- as.data.frame(multiScaling2())
      multiscaling$color <- as.vector(selectedDynamicColor2())
      my_colors <- as.vector(unique(multiscaling$color))
      my_colors <- my_colors[order(my_colors)]
      svg("Multiscaling_1_2.svg", width = 10, height = 10)
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale1, y = scale2, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 1",
                   y = "Scaling Dimension 2")
            + scale_color_manual(values=my_colors)) 
      dev.off()
      
      fs <- c(fs, "Multiscaling_1_2.svg")
      
      svg("Multiscaling_2_3.svg", width = 10, height = 10)
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale2, y = scale3, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 2",
                   y = "Scaling Dimension 3")
            + scale_color_manual(values=my_colors))
      dev.off()
      
      fs <- c(fs, "Multiscaling_2_3.svg")
      
      svg("Multiscaling_1_3.svg", width = 10, height = 10)
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale1, y = scale3, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 1",
                   y = "Scaling Dimension 3")
            + scale_color_manual(values=my_colors))
      dev.off()
      
      fs <- c(fs, "Multiscaling_1_3.svg")
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

######################

output$downloadEigenD3 <- downloadHandler(
  filename = function(){
    "Modules_Eigenvalue.csv"
  },
  content = function (filename){
    write.csv(selectedMEs3(), filename)
  }
)

output$downloadModuleD3 <- downloadHandler(
  filename = function(){
    "Modules_Size.csv"
  },
  content = function (filename){
    write.csv(as.data.frame(table(selectedDynamicColor3())), filename)
  }
)

output$Download_Network_Creation_dataset3 <- downloadHandler(
  filename = function(){
    paste("NETWORK_CREATION_dataset3.zip")
  },
  content = function (filename){
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    if (input$pdf_or_svg_p3_dataset3 == "pdf"){
      pdf("Soft_Power.pdf", width = 10, height = 10)
      print(ggplot(data =sft3()$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq )) + 
              geom_text(aes(label = Power, color = "tomato3"), show.legend =  FALSE) +
              geom_hline(yintercept = 0.9) +
              labs(y = "Scale Free Topology", x = "Soft Threshold (Power)", title = "Scale Independence") +
              theme(legend.position = "none"))
      dev.off()
      
      fs <- c(fs, "Soft_Power.pdf")
      pdf("Module_dendrogramme.pdf", width = 10, height = 10)
      print(ggdendrogram(selectedMETree3()) + 
              labs(y = "Height", title = "Clustering of module eigengenes"))
      dev.off()
      pdf("Multiscaling.pdf", width = 10, height = 10)
      multiscaling <- as.data.frame(multiScaling3())
      multiscaling$color <- as.vector(selectedDynamicColor3())
      my_colors <- as.vector(unique(multiscaling$color))
      my_colors <- my_colors[order(my_colors)]
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale1, y = scale2, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 1",
                   y = "Scaling Dimension 2")
            + scale_color_manual(values=my_colors)) 
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale2, y = scale3, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 2",
                   y = "Scaling Dimension 3")
            + scale_color_manual(values=my_colors))
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale1, y = scale3, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 1",
                   y = "Scaling Dimension 3")
            + scale_color_manual(values=my_colors))
      dev.off()
      
      fs <- c(fs, "Multiscaling.pdf")
      fs <- c(fs, "Module_dendrogramme.pdf")
      pdf("gene_dendrogramme_colors.pdf", width = 10, height = 10)
      print(plotDendroAndColors(selectedTree3(), selectedDynamicColor3(), "Dynamic Tree Cut",
                                dendroLabels = FALSE, hang = 0.03,
                                addGuide = TRUE, guideHang = 0.05,
                                main = "Gene dendrogram and module colors"))
      dev.off()
      fs <- c(fs, "gene_dendrogramme_colors.pdf")
    }else{
      svg("Soft_Power.svg", width = 10, height = 10)
      print(ggplot(data =sft3()$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq )) + 
              geom_text(aes(label = Power, color = "tomato3"), show.legend =  FALSE) +
              geom_hline(yintercept = 0.9) +
              labs(y = "Scale Free Topology", x = "Soft Threshold (Power)", title = "Scale Independence") +
              theme(legend.position = "none"))
      dev.off()
      
      fs <- c(fs, "Soft_Power.svg")
      svg("Module_dendrogramme.svg", width = 10, height = 10)
      print(ggdendrogram(selectedMETree3()) + 
              labs(y = "Height", title = "Clustering of module eigengenes"))
      dev.off()
      multiscaling <- as.data.frame(multiScaling3())
      multiscaling$color <- as.vector(selectedDynamicColor3())
      my_colors <- as.vector(unique(multiscaling$color))
      my_colors <- my_colors[order(my_colors)]
      svg("Multiscaling_1_2.svg", width = 10, height = 10)
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale1, y = scale2, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 1",
                   y = "Scaling Dimension 2")
            + scale_color_manual(values=my_colors)) 
      dev.off()
      
      fs <- c(fs, "Multiscaling_1_2.svg")
      
      svg("Multiscaling_2_3.svg", width = 10, height = 10)
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale2, y = scale3, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 2",
                   y = "Scaling Dimension 3")
            + scale_color_manual(values=my_colors))
      dev.off()
      
      fs <- c(fs, "Multiscaling_2_3.svg")
      
      svg("Multiscaling_1_3.svg", width = 10, height = 10)
      print(ggplot(data=as.data.frame(multiscaling)) 
            + geom_point(aes(x = scale1, y = scale3, color = color))
            + labs(title = "MDS plot", x = "Scaling Dimension 1",
                   y = "Scaling Dimension 3")
            + scale_color_manual(values=my_colors))
      dev.off()
      
      fs <- c(fs, "Multiscaling_1_3.svg")
      fs <- c(fs, "Module_dendrogramme.svg")
      svg("gene_dendrogramme_colors.svg", width = 10, height = 10)
      print(plotDendroAndColors(selectedTree3(), selectedDynamicColor3(), "Dynamic Tree Cut",
                                dendroLabels = FALSE, hang = 0.03,
                                addGuide = TRUE, guideHang = 0.05,
                                main = "Gene dendrogram and module colors"))
      dev.off()
      fs <- c(fs, "gene_dendrogramme_colors.svg")
    }
    
    zip(zipfile=filename, files=fs)
    
  },
  contentType = "application/zip")


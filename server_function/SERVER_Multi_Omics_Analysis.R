#################################################
#### SERVER FUNCTIONS : MULTI-OMICS ANALYSIS ####
#################################################
# This script gathers the functions related to the multi-omics analysis tabulation of MiBiOmics. It contains all the
# functions needed for the co-inertia analysis, procrustes analysis and the association of modules extracted from
# WGCNA.

source("function.R")
#### REACTIVE OBJECTS ####

mcoia <- reactive({
  if (input$Omic3){
    selectD1 <- selected_moduleOTU()
    selectD2 <- selectedMetab()
    selectD3 <- selectedD3()
    df1_df2_df3 <- list()
    df1_df2_df3[["df1"]] <- t(as.data.frame(selectD1))
    df1_df2_df3[["df2"]] <- t(as.data.frame(selectD2))
    df1_df2_df3[["df3"]] <- t(as.data.frame(selectD3))
    rownames(df1_df2_df3[["df1"]]) <- paste("a_", rownames(df1_df2_df3[["df1"]]), sep = "")
    rownames(df1_df2_df3[["df2"]]) <- paste("b_", rownames(df1_df2_df3[["df2"]]), sep = "")
    rownames(df1_df2_df3[["df3"]]) <- paste("c_", rownames(df1_df2_df3[["df3"]]), sep = "")

    mcia(df1_df2_df3, cia.scan = TRUE, cia.nf = 2)
  }else{
    selectMetab <- selectedMetab()
    df1_df2 <- list()
    df1_df2[["df1"]] <- t(as.data.frame(selected_moduleOTU()))
    df1_df2[["df2"]] <- t(as.data.frame(selectMetab))
    rownames(df1_df2[["df1"]]) <- paste("a_", rownames(df1_df2[["df1"]]), sep = "")
    rownames(df1_df2[["df2"]]) <- paste("b_", rownames(df1_df2[["df2"]]), sep = "")
    mcia(df1_df2, nsc = F)
  }
})


selected_sampleInfo <- reactive({
  conditionCol <- input$SelectVariable2
  if (input$Omic3){
    vector_info1 <- c(rep(sampleAnnot_3()[,conditionCol], 3))
  }else{
    vector_info1 <- c(rep(sampleAnnot_3()[,conditionCol], 2))
  }
  sample_info_mcoia <- data.frame(as.data.frame(mcoia()[["mcoa"]]$Tl1), as.data.frame(vector_info1))
  sample_info_mcoia
})

LegendDF <- reactive({
  sample_info_mcoia <- selected_sampleInfo()
  sample_info_mcoia$DF <- substr(rownames(sample_info_mcoia), nchar(rownames(sample_info_mcoia))-2, nchar(rownames(sample_info_mcoia)))

  if (input$Omic3){
    PlotLegend <- ggplot(data = sample_info_mcoia) +
      geom_point(aes_string(x=colnames(sample_info_mcoia)[1], y=colnames(sample_info_mcoia)[2], shape = "DF"))+
      scale_shape_manual(name="Omics", labels=c(input$CountingT1, input$OmicTable, input$OmicTable3), values=c(15, 16, 17))
  }else{
    PlotLegend <- ggplot(data = sample_info_mcoia) +
      geom_point(aes_string(x=colnames(sample_info_mcoia)[1], y=colnames(sample_info_mcoia)[2], shape = "DF"))+
      scale_shape_manual(name="Omics", labels=c(input$CountingT1, input$OmicTable), values=c(15, 16))
  }

  myLegend <- get_legend(PlotLegend)
  myLegend
})

selected_sampleInfo_df1 <- reactive({
  if (input$Omic3){
    sample_info_mcoia_df1 <-selected_sampleInfo()[1:(nrow(selected_sampleInfo())/3),1:2]
  }else{
    sample_info_mcoia_df1 <-selected_sampleInfo()[1:(nrow(selected_sampleInfo())/2),1:2]
  }
  sample_info_mcoia_df1
})


selected_sampleInfo_df2 <- reactive({
  if (input$Omic3){
    sample_info_mcoia_df2 <-selected_sampleInfo()[(nrow(selected_sampleInfo())/3+1):(nrow(selected_sampleInfo())*(2/3)),1:2]
  }else{
    sample_info_mcoia_df2 <-selected_sampleInfo()[(nrow(selected_sampleInfo())/2+1):(nrow(selected_sampleInfo())),1:2]
  }
  sample_info_mcoia_df2
})

selected_sampleInfo_df3 <- reactive({
  if (input$Omic3){
    sample_info_mcoia_df3 <-selected_sampleInfo()[(nrow(selected_sampleInfo())*(2/3)+1):(nrow(selected_sampleInfo())),1:2]
  }else{
    sample_info_mcoia_df3 <- 0
  }
  sample_info_mcoia_df3
})


selected_sampleInfo_all <- reactive({
  if (input$Omic3){
    sample_info_mcoia_all <- data.frame(selected_sampleInfo_df1(), selected_sampleInfo_df2(), selected_sampleInfo_df3(), as.data.frame(selected_sampleInfo()))
    colnames(sample_info_mcoia_all) <- c("col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8", "col9")
  }else{
    sample_info_mcoia_all <- data.frame(selected_sampleInfo_df1(), selected_sampleInfo_df2(), selected_sampleInfo())
    colnames(sample_info_mcoia_all) <- c("col1", "col2", "col3", "col4", "col5", "col6", "col7")
  }
  sample_info_mcoia_all
})


### A dataframe containing the information about the length of each vector in the co-inertia
distance_mcoia_df <- reactive({
  validate(
    need(!input$Omic3, 'Only for multi-omics analysis with two datasets')
  )
  distance_coia.v <- c()
  print(selected_sampleInfo_all())
  for (i in 1:nrow(selected_sampleInfo_all())){
    # Calculate distance between two point in co-inertia
    dist_coia <- sqrt((selected_sampleInfo_all()[i, 4] - selected_sampleInfo_all()[i, 2])^2 + (selected_sampleInfo_all()[i, 3] - selected_sampleInfo_all()[i, 1])^2)
    distance_coia.v <- c(distance_coia.v, dist_coia)
  }
  distance_mcoia.df <- data.frame(distance_coia.v)
  rownames(distance_mcoia.df) <- rownames(sampleAnnot_3())
  distance_mcoia.df$sampleName <- rownames(sampleAnnot_3())
  distance_mcoia.df

})

distTree_coia <- reactive({
  validate(
    need(!input$Omic3, 'Only for multi-omics analysis with two datasets')
  )
  hclust(dist(distance_mcoia_df()), method = "average")
})

distance_mcoia_df_supp_info <- reactive({
  validate(
    need(!input$Omic3, 'Only for multi-omics analysis with two datasets')
  )
  distance_mcoia.df <- distance_mcoia_df()
  distance_mcoia.df$color <- sampleAnnot_3()[,input$SelectVariable2]

  distance_mcoia.df <- distance_mcoia.df[order(distance_mcoia.df$distance_coia.v),]
  distance_mcoia.df$sampleName <- as.character(distance_mcoia.df$sampleName)
  distance_mcoia.df$sampleName <- factor(distance_mcoia.df$sampleName, levels= unique(distance_mcoia.df$sampleName))
  distance_mcoia.df
})

bivariate_plot <- reactive({
  validate(
    need(!input$Omic3, 'Only for multi-omics analysis with two datasets')
  )
  bivariate_plot <- data.frame(mcoia()[["mcoa"]]$Tl1[1:(nrow(selected_sampleInfo())/2), 1], mcoia()[["mcoa"]]$Tl1[(nrow(selected_sampleInfo())/2+1):(nrow(selected_sampleInfo())), 1])

  colnames(bivariate_plot) <- c(paste("Coinertia_", input$CountingT1, "_axis", sep =""), paste("Coinertia_", input$OmicTable, "D2_axis", sep=""))

  bivariate_plot$sampleName <- rownames(sampleAnnot_3())
  bivariate_plot
})

RV_score <- reactive({
  mcoia()[["mcoa"]]$RV[1,2]
})

selected_axis1_drivers <- reactive({
  axis1_drivers <- data.frame(Variables_PosX = character(), Variables_PosY = character())
  mco <- mcoia()
  for (i in 1:nrow(mco[["mcoa"]]$Tco)){
    if (input$ShowDrivers_axis == "axis1"){
      if (mco[["mcoa"]]$Tco$SV1[i] < quantile(mco[["mcoa"]]$Tco$SV1, 0.01)){
        axis1_drivers <- rbind(axis1_drivers, mco[["mcoa"]]$Tco[i,])
      }
      if (mco[["mcoa"]]$Tco$SV1[i] > quantile(mco[["mcoa"]]$Tco$SV1, 0.99)){
        axis1_drivers <- rbind(axis1_drivers, mco[["mcoa"]]$Tco[i,])
      }
    }else{
      if (mco[["mcoa"]]$Tco$SV2[i] < quantile(mco[["mcoa"]]$Tco$SV2, 0.01)){
        axis1_drivers <- rbind(axis1_drivers, mco[["mcoa"]]$Tco[i,])
      }
      if (mco[["mcoa"]]$Tco$SV2[i] > quantile(mco[["mcoa"]]$Tco$SV2, 0.99)){
        axis1_drivers <- rbind(axis1_drivers, mco[["mcoa"]]$Tco[i,])
      }
    }
  }
  axis_drivers <- axis1_drivers
  if (input$Omic3){
    name_axis_drivers <- rename_axis_drivers_simplified(axis1_drivers)
    feature_type_axis_1 <- c()
    spec_axis_1 <- c()
    the_name <- c()
    drivers <- data.frame("old" = rownames(axis1_drivers), "new" = name_axis_drivers)
    for (i in 1:nrow(drivers)){
      if (substr(drivers$old[i], 1, 1) == "a"){
        feature_type_axis_1 <- c(feature_type_axis_1, input$CountingT1)
        if (input$CountingT1 == "OTUs" && input$TaxonFile1){
          spec_axis_1 <- c(spec_axis_1, as.character(taxTable()[as.character(drivers$new[i]), input$selectTaxonomy]))
        }else{
          spec_axis_1 <- c(spec_axis_1, as.character(drivers$new[i]))
        }

      }else{
        if (substr(drivers$old[i], 1, 1) == "b"){
          feature_type_axis_1 <- c(feature_type_axis_1, input$OmicTable)
          if (input$OmicTable == "OTUs" && input$TaxonFile1){
            spec_axis_1 <- c(spec_axis_1, as.character(taxTable()[as.character(drivers$new[i]), input$selectTaxonomy]))
          }else{
            spec_axis_1 <- c(spec_axis_1, as.character(drivers$new[i]))
          }

        }else{
          if (substr(drivers$old[i], 1, 1) == "c"){
            feature_type_axis_1 <- c(feature_type_axis_1, input$OmicTable3)
            if (input$OmicTable3 == "OTUs" && input$TaxonFile1){
              spec_axis_1 <- c(spec_axis_1, as.character(taxTable()[as.character(drivers$new[i]), input$selectTaxonomy]))
            }else{
              spec_axis_1 <- c(spec_axis_1, as.character(drivers$new[i]))
            }

          }else{
            feature_type_axis_1 <- c(feature_type_axis_1, "unkwown")
            spec_axis_1 <- c(spec_axis_1, as.character(drivers$new[i]))
          }
        }
      }
    }
  }else{
    name_axis_drivers <- rename_axis_drivers_simplified(axis1_drivers)
    feature_type_axis_1 <- c()
    spec_axis_1 <- c()
    the_name <- c()
    drivers <- data.frame("old" = rownames(axis1_drivers), "new" = name_axis_drivers)
    for (i in 1:nrow(drivers)){
      if (substr(drivers$old[i], 1, 1) == "a"){

        feature_type_axis_1 <- c(feature_type_axis_1, input$CountingT1)
        if (input$CountingT1 == "OTUs" && input$TaxonFile1){
          spec_axis_1 <- c(spec_axis_1, as.character(taxTable()[as.character(drivers$new[i]), input$selectTaxonomy]))
        }else{
          spec_axis_1 <- c(spec_axis_1, as.character(drivers$new[i]))
        }



      }else{
        if (substr(drivers$old[i], 1, 1) == "b"){

          feature_type_axis_1 <- c(feature_type_axis_1, input$OmicTable)
          if (input$OmicTable == "OTUs" && input$TaxonFile1){
            spec_axis_1 <- c(spec_axis_1, as.character(taxTable()[as.character(drivers$new[i]), input$selectTaxonomy]))
          }else{
            spec_axis_1 <- c(spec_axis_1, as.character(drivers$new[i]))
          }

        }else{
          feature_type_axis_1 <- c(feature_type_axis_1, "unkwown")
          spec_axis_1 <- c(spec_axis_1, as.character(drivers$new[i]))
        }
      }
    }
  }

  axis1_drivers$feature <- feature_type_axis_1
  axis1_drivers$spec <- spec_axis_1
  axis1_drivers$name <- the_name
  axis1_drivers
})


#### PROCRUSTES


dudi_exprDat <- reactive({

  # if (input$selectModule1_p5 == "General"){
  dudi_exprDat <- dudi.pca(exprDat_3(), scannf = FALSE, nf = 2 )
  # }else{
  #   Module_OTU <- exprDat_3()[,modGenes]
  #   dudi_exprDat <- dudi.pca(Module_OTU, scannf = FALSE, nf = 2 )
  # }

  dudi_exprDat
})


dudi_exprDatSec <- reactive({

  dudi_exprDatSec <- dudi.pca(selectedMetab(), scannf = FALSE, nf = 2 )

  dudi_exprDatSec
})



PA_object <- reactive({
  validate(
    need(!input$Omic3, 'Only for multi-omics analysis with two datasets')
  )
  PA_object <- procrustes(dudi_exprDat(), dudi_exprDatSec(), scale = TRUE)
  PA_object
})


sample_info_PA <- reactive({
  validate(
    need(!input$Omic3, 'Only for multi-omics analysis with two datasets')
  )
  sample_info_X_PA <- data.frame(PA_object()[["X"]], sampleAnnot_3())
  colnames(sample_info_X_PA) <- c("Axis1", "Axis2", colnames(sampleAnnot_3()))
  sample_info_Y_PA <- data.frame(PA_object()[["Yrot"]], sampleAnnot_3())
  colnames(sample_info_Y_PA) <- c("Axis.1", "Axis.2", colnames(sampleAnnot_3()))
  sample_info_PA <- data.frame(sample_info_X_PA, sample_info_Y_PA)
  sample_info_PA
})


#### PRE - HEATMAP ####



corr_MEs_D1D2 <- reactive({
  MEs <- selectedMEs()
  MEs_2 <- selectedMEs2()
  MEs <- MEs[which(rownames(MEs) %in% rownames(MEs_2)),]
  MEs_2 <- MEs_2[which(rownames(MEs_2) %in% rownames(MEs)),]
  corr_MEs_D1D2 <- cor(MEs, MEs_2, use = "p", method = "spearman")
  corr_MEs_D1D2 <- as.data.frame(corr_MEs_D1D2)
  row.order <- hclust(dist(corr_MEs_D1D2, method = "euclidean"), method = "ward.D")$order
  col.order <- hclust(dist(t(corr_MEs_D1D2), method = "euclidean"), method = "ward.D")$order
  corr_MEs_D1D2 <- corr_MEs_D1D2[row.order, col.order]
  t(corr_MEs_D1D2)
})

corr_MEs_D2D3 <- reactive({
  MEs <- selectedMEs2()
  MEs_2 <- selectedMEs3()
  MEs <- MEs[which(rownames(MEs) %in% rownames(MEs_2)),]
  MEs_2 <- MEs_2[which(rownames(MEs_2) %in% rownames(MEs)),]
  corr_MEs_D1D2 <- cor(MEs, MEs_2, use = "p", method = "spearman")
  corr_MEs_D1D2 <- as.data.frame(corr_MEs_D1D2)
  row.order <- hclust(dist(corr_MEs_D1D2, method = "euclidean"), method = "ward.D")$order
  col.order <- hclust(dist(t(corr_MEs_D1D2), method = "euclidean"), method = "ward.D")$order
  corr_MEs_D1D2 <- corr_MEs_D1D2[row.order, col.order]
  t(corr_MEs_D1D2)
})

corr_MEs_D1D3 <- reactive({
  MEs <- selectedMEs()
  MEs_2 <- selectedMEs3()
  MEs <- MEs[which(rownames(MEs) %in% rownames(MEs_2)),]
  MEs_2 <- MEs_2[which(rownames(MEs_2) %in% rownames(MEs)),]
  corr_MEs_D1D2 <- cor(MEs, MEs_2, use = "p", method = "spearman")
  corr_MEs_D1D2 <- as.data.frame(corr_MEs_D1D2)
  row.order <- hclust(dist(corr_MEs_D1D2, method = "euclidean"), method = "ward.D")$order
  col.order <- hclust(dist(t(corr_MEs_D1D2), method = "euclidean"), method = "ward.D")$order
  corr_MEs_D1D2 <- corr_MEs_D1D2[row.order, col.order]
  t(corr_MEs_D1D2)
})


#### HIVE PLOT ####

hive_data <- reactive({
  myCorr_D1_D2 <- list()

  myCorr_D1_D2[[1]] <- corr_MEs_D1D2()
  myCorr_D1_D2[[2]] <- p_val_MEs_D1D2()

  rownames(myCorr_D1_D2[[1]]) <- paste(substr(rownames(myCorr_D1_D2[[1]]), 1, nchar(rownames(myCorr_D1_D2[[1]]))), "_DF2" ,sep = "")
  rownames(myCorr_D1_D2[[2]]) <- paste(substr(rownames(myCorr_D1_D2[[2]]), 1, nchar(rownames(myCorr_D1_D2[[2]]))), "_DF2", sep ="")
  colnames(myCorr_D1_D2[[1]]) <- paste(substr(colnames(myCorr_D1_D2[[1]]), 1, nchar(colnames(myCorr_D1_D2[[1]]))), "_DF1" ,sep = "")
  colnames(myCorr_D1_D2[[2]]) <- paste(substr(colnames(myCorr_D1_D2[[2]]), 1, nchar(colnames(myCorr_D1_D2[[2]]))), "_DF1", sep ="")

  WGCNA_1 <- list()
  WGCNA_1[[5]] <- selectedMEs()
  WGCNA_2 <- list()
  WGCNA_2[[5]] <- selectedMEs2()

  my_Corrs <- list()
  my_Corrs[[1]] <- myCorr_D1_D2
  my_Annot <- list()
  my_Annot[[1]] <- sampleAnnot_WGCNA()
  my_Annot[[2]] <- sampleAnnot_sec_WGCNA()

  if (input$Omic3){
    myCorr_D1_D3 <- list()
    myCorr_D2_D3 <- list()
    myCorr_D1_D3[[1]] <- corr_MEs_D1D3()
    myCorr_D1_D3[[2]] <- p_val_MEs_D1D3()
    myCorr_D2_D3[[1]] <- corr_MEs_D2D3()
    myCorr_D2_D3[[2]] <- p_val_MEs_D2D3()
    rownames(myCorr_D1_D3[[1]]) <- paste(substr(rownames(myCorr_D1_D3[[1]]), 1, nchar(rownames(myCorr_D1_D3[[1]]))), "_DF3" ,sep = "")
    rownames(myCorr_D1_D3[[2]]) <- paste(substr(rownames(myCorr_D1_D3[[2]]), 1, nchar(rownames(myCorr_D1_D3[[2]]))), "_DF3", sep ="")
    colnames(myCorr_D1_D3[[1]]) <- paste(substr(colnames(myCorr_D1_D3[[1]]), 1, nchar(colnames(myCorr_D1_D3[[1]]))), "_DF1" ,sep = "")
    colnames(myCorr_D1_D3[[2]]) <- paste(substr(colnames(myCorr_D1_D3[[2]]), 1, nchar(colnames(myCorr_D1_D3[[2]]))), "_DF1", sep ="")

    colnames(myCorr_D2_D3[[1]]) <- paste(substr(colnames(myCorr_D2_D3[[1]]), 1, nchar(colnames(myCorr_D2_D3[[1]]))), "_DF2" ,sep = "")
    colnames(myCorr_D2_D3[[2]]) <- paste(substr(colnames(myCorr_D2_D3[[2]]), 1, nchar(colnames(myCorr_D2_D3[[2]]))), "_DF2", sep ="")
    rownames(myCorr_D2_D3[[1]]) <- paste(substr(rownames(myCorr_D2_D3[[1]]), 1, nchar(rownames(myCorr_D2_D3[[1]]))), "_DF3" ,sep = "")
    rownames(myCorr_D2_D3[[2]]) <- paste(substr(rownames(myCorr_D2_D3[[2]]), 1, nchar(rownames(myCorr_D2_D3[[2]]))), "_DF3", sep ="")

    WGCNA_3 <- list()
    WGCNA_3[[5]] <- selectedMEs3()

    my_Corrs[[2]] <- myCorr_D1_D3
    my_Corrs[[3]] <- myCorr_D2_D3
    my_Annot[[3]] <- sampleAnnot_ter_WGCNA()
    my_hive <- hive_myLayers(WGCNA_1, WGCNA_2, WGCNA_3, myCorrs = my_Corrs, myAnnots = my_Annot, trait = input$traitHive)
    print(my_hive)
    print(nrow(edges))
    validate(
      need(nrow(my_hive$edges) != 0, "No significant interaction between modules and external trait")
    )
    validate(
      need(my_hive$edges != 0, "No significant interaction between modules and external trait")
    )
  }else{
    my_hive <- hive_my2Layers(WGCNA_1, WGCNA_2, myCorr = myCorr_D1_D2, myAnnots = my_Annot, trait = input$traitHive)
    print(my_hive)

    validate(
      need(my_hive != 0, "No significant interaction between modules and external trait")
    )
    validate(
      need(nrow(my_hive$edges) != 0, "No significant interaction between modules and external trait")
    )
  }
  my_hive

})

# hive_data_keystone <- reactive({
#   myCorr_D1_D2 <- list()
#
#   myCorr_D1_D2[[1]] <- corr_MEs_D1D2()
#   myCorr_D1_D2[[2]] <- p_val_MEs_D1D2()
#
#   rownames(myCorr_D1_D2[[1]]) <- paste(substr(rownames(myCorr_D1_D2[[1]]), 1, nchar(rownames(myCorr_D1_D2[[1]]))), "_DF2" ,sep = "")
#   rownames(myCorr_D1_D2[[2]]) <- paste(substr(rownames(myCorr_D1_D2[[2]]), 1, nchar(rownames(myCorr_D1_D2[[2]]))), "_DF2", sep ="")
#   colnames(myCorr_D1_D2[[1]]) <- paste(substr(colnames(myCorr_D1_D2[[1]]), 1, nchar(colnames(myCorr_D1_D2[[1]]))), "_DF1" ,sep = "")
#   colnames(myCorr_D1_D2[[2]]) <- paste(substr(colnames(myCorr_D1_D2[[2]]), 1, nchar(colnames(myCorr_D1_D2[[2]]))), "_DF1", sep ="")
#
#   WGCNA_1 <- list()
#   WGCNA_1[[5]] <- selectedMEs()
#   WGCNA_2 <- list()
#   WGCNA_2[[5]] <- selectedMEs2()
#
#   my_Corrs <- list()
#   my_Corrs[[1]] <- myCorr_D1_D2
#   my_Annot <- list()
#   my_Annot[[1]] <- sampleAnnot_WGCNA()
#   my_Annot[[2]] <- sampleAnnot_sec_WGCNA()
#
#   if (input$Omic3){
#     myCorr_D1_D3 <- list()
#     myCorr_D2_D3 <- list()
#     myCorr_D1_D3[[1]] <- corr_MEs_D1D3()
#     myCorr_D1_D3[[2]] <- p_val_MEs_D1D3()
#     myCorr_D2_D3[[1]] <- corr_MEs_D2D3()
#     myCorr_D2_D3[[2]] <- p_val_MEs_D2D3()
#     rownames(myCorr_D1_D3[[1]]) <- paste(substr(rownames(myCorr_D1_D3[[1]]), 1, nchar(rownames(myCorr_D1_D3[[1]]))), "_DF3" ,sep = "")
#     rownames(myCorr_D1_D3[[2]]) <- paste(substr(rownames(myCorr_D1_D3[[2]]), 1, nchar(rownames(myCorr_D1_D3[[2]]))), "_DF3", sep ="")
#     colnames(myCorr_D1_D3[[1]]) <- paste(substr(colnames(myCorr_D1_D3[[1]]), 1, nchar(colnames(myCorr_D1_D3[[1]]))), "_DF1" ,sep = "")
#     colnames(myCorr_D1_D3[[2]]) <- paste(substr(colnames(myCorr_D1_D3[[2]]), 1, nchar(colnames(myCorr_D1_D3[[2]]))), "_DF1", sep ="")
#
#     colnames(myCorr_D2_D3[[1]]) <- paste(substr(colnames(myCorr_D2_D3[[1]]), 1, nchar(colnames(myCorr_D2_D3[[1]]))), "_DF2" ,sep = "")
#     colnames(myCorr_D2_D3[[2]]) <- paste(substr(colnames(myCorr_D2_D3[[2]]), 1, nchar(colnames(myCorr_D2_D3[[2]]))), "_DF2", sep ="")
#     rownames(myCorr_D2_D3[[1]]) <- paste(substr(rownames(myCorr_D2_D3[[1]]), 1, nchar(rownames(myCorr_D2_D3[[1]]))), "_DF3" ,sep = "")
#     rownames(myCorr_D2_D3[[2]]) <- paste(substr(rownames(myCorr_D2_D3[[2]]), 1, nchar(rownames(myCorr_D2_D3[[2]]))), "_DF3", sep ="")
#
#     WGCNA_3 <- list()
#     WGCNA_3[[5]] <- selectedMEs3()
#
#     my_Corrs[[2]] <- myCorr_D1_D3
#     my_Corrs[[3]] <- myCorr_D2_D3
#     my_Annot[[3]] <- sampleAnnot_ter_WGCNA()
#     my_hive <- hive_myLayers(WGCNA_1, WGCNA_2, WGCNA_3, myCorrs = my_Corrs, myAnnots = my_Annot, trait = input$traitHive, networkOrdered = TRUE)
#     print(my_hive)
#     print(nrow(edges))
#     validate(
#       need(nrow(my_hive$edges) != 0, "No significant interaction between modules and external trait")
#     )
#     validate(
#       need(my_hive$edges != 0, "No significant interaction between modules and external trait")
#     )
#   }else{
#     my_hive <- hive_my2Layers(WGCNA_1, WGCNA_2, myCorr = myCorr_D1_D2, myAnnots = my_Annot, trait = input$traitHive, networkOrdered = TRUE)
#     print(my_hive)
#
#     validate(
#       need(my_hive != 0, "No significant interaction between modules and external trait")
#     )
#     validate(
#       need(nrow(my_hive$edges) != 0, "No significant interaction between modules and external trait")
#     )
#   }
#   my_hive
#
# })



##########################
p_val_MEs_D1D2 <- reactive({
  MEs <- selectedMEs()
  MEs2 <- selectedMEs2()
  MEs <- MEs[which(rownames(MEs) %in% rownames(MEs2)),]
  MEs2 <- MEs2[which(rownames(MEs2) %in% rownames(MEs)),]
  pvalue_MEs_D1_D2 <- corPvalueStudent(cor(MEs, MEs2, use = "p", method = "spearman"), nrow(MEs))
  pvalue_MEs_D1_D2 <- as.data.frame(pvalue_MEs_D1_D2)
  row.order <- hclust(dist(pvalue_MEs_D1_D2, method = "euclidean"), method = "ward.D")$order
  col.order <- hclust(dist(t(pvalue_MEs_D1_D2), method = "euclidean"), method = "ward.D")$order
  pvalue_MEs_D1_D2 <- pvalue_MEs_D1_D2[row.order, col.order]
  t(pvalue_MEs_D1_D2)
})

p_val_MEs_D2D3 <- reactive({
  MEs3 <- selectedMEs3()
  MEs2 <- selectedMEs2()
  MEs3 <- MEs3[which(rownames(MEs3) %in% rownames(MEs2)),]
  MEs2 <- MEs2[which(rownames(MEs2) %in% rownames(MEs3)),]
  pvalue_MEs_D1_D2 <- corPvalueStudent(cor(MEs2, MEs3, use = "p", method = "spearman"), nrow(MEs2))
  pvalue_MEs_D1_D2 <- as.data.frame(pvalue_MEs_D1_D2)
  row.order <- hclust(dist(pvalue_MEs_D1_D2, method = "euclidean"), method = "ward.D")$order
  col.order <- hclust(dist(t(pvalue_MEs_D1_D2), method = "euclidean"), method = "ward.D")$order
  pvalue_MEs_D1_D2 <- pvalue_MEs_D1_D2[row.order, col.order]
  t(pvalue_MEs_D1_D2)
})

p_val_MEs_D1D3 <- reactive({
  MEs3 <- selectedMEs3()
  MEs <- selectedMEs()
  MEs3 <- MEs3[which(rownames(MEs3) %in% rownames(MEs)),]
  MEs <- MEs[which(rownames(MEs) %in% rownames(MEs3)),]
  pvalue_MEs_D1_D2 <- corPvalueStudent(cor(MEs, MEs3, use = "p", method = "spearman"), nrow(MEs))
  pvalue_MEs_D1_D2 <- as.data.frame(pvalue_MEs_D1_D2)
  row.order <- hclust(dist(pvalue_MEs_D1_D2, method = "euclidean"), method = "ward.D")$order
  col.order <- hclust(dist(t(pvalue_MEs_D1_D2), method = "euclidean"), method = "ward.D")$order
  pvalue_MEs_D1_D2 <- pvalue_MEs_D1_D2[row.order, col.order]
  t(pvalue_MEs_D1_D2)
})

corr_expr <- reactive({
  # print(exprDatSec_51())
  # validate(
  #   need(any(rowSums(exprDat_41()) == 0), "Modules' variables are not expressed in your samples. Please transform your data to avoid zeros.")
  # )
  # validate(
  #   need(any(rowSums(exprDatSec_51()) == 0), "Modules' variables are not expressed in your samples. Please transform your data to avoid zeros.")
  # )

  corr_expr_1 <- cor(exprDat_41(), exprDatSec_51(), use = "p", method = "pearson")
  corr_expr_1 <- as.data.frame(corr_expr_1)
  rownames(corr_expr_1) <- paste("DF1_", rownames(corr_expr_1), sep = "")
  colnames(corr_expr_1) <- paste("DF2_", colnames(corr_expr_1), sep = "")
  corr_expr_1[is.na(corr_expr_1)] <- 0
  corr_expr_1 <- corr_expr_1[, which(colSums(corr_expr_1) != 0)]
  row.order <- hclust(dist(corr_expr_1, method = "euclidean"), method = "ward.D")$order
  col.order <- hclust(dist(t(corr_expr_1), method = "euclidean"), method = "ward.D")$order
  corr_expr_new <- corr_expr_1[row.order, col.order]
  t(corr_expr_new)
})

corr_expr23 <- reactive({
  # validate(
  #   need(any(rowSums(exprDatSec_52()) == 0), "Your second omics dataset must be filtrated to remove low counts")
  # )
  # validate(
  #   need(any(rowSums(exprDatTer_52()) == 0), "Your third omics dataset must be filtrated to remove low counts")
  # )
  corr_expr_1 <- cor(exprDatSec_52(), exprDatTer_52(), use = "p", method = "pearson")
  corr_expr_1 <- as.data.frame(corr_expr_1)
  rownames(corr_expr_1) <- paste("DF2_", rownames(corr_expr_1), sep = "")
  colnames(corr_expr_1) <- paste("DF3_", colnames(corr_expr_1), sep = "")
  corr_expr_1[is.na(corr_expr_1)] <- 0
  corr_expr_1 <- corr_expr_1[, which(colSums(corr_expr_1) != 0)]
  row.order <- hclust(dist(corr_expr_1, method = "euclidean"), method = "ward.D")$order
  col.order <- hclust(dist(t(corr_expr_1), method = "euclidean"), method = "ward.D")$order
  corr_expr_new <- corr_expr_1[row.order, col.order]
  t(corr_expr_new)
})

corr_expr13 <- reactive({
  # validate(
  #   need(any(rowSums(exprDat_43()) == 0), "Your first omics dataset must be filtrated to remove low counts")
  # )
  # validate(
  #   need(any(rowSums(exprDatTer_53()) == 0), "Your third omics dataset must be filtrated to remove low counts")
  # )
  print(exprDat_43())
  print(exprDatTer_53())
  corr_expr_1 <- cor(exprDat_43(), exprDatTer_53(), use = "p", method = "pearson")
  corr_expr_1 <- as.data.frame(corr_expr_1)
  rownames(corr_expr_1) <- paste("DF1_", rownames(corr_expr_1), sep = "")
  colnames(corr_expr_1) <- paste("DF3_", colnames(corr_expr_1), sep = "")
  corr_expr_1[is.na(corr_expr_1)] <- 0
  corr_expr_1 <- corr_expr_1[, which(colSums(corr_expr_1) != 0)]
  row.order <- hclust(dist(corr_expr_1, method = "euclidean"), method = "ward.D")$order
  col.order <- hclust(dist(t(corr_expr_1), method = "euclidean"), method = "ward.D")$order
  corr_expr_new <- corr_expr_1[row.order, col.order]
  t(corr_expr_new)
})


family_group12_1 <- reactive({
  if (input$CountingT1 == "OTUs"){
    taxon_group <- c()
    for (i in 1:ncol(corr_expr())){
      taxon_group <- c(taxon_group, (taxTable2()[which(taxTable2()[["rn"]] == colnames(corr_expr())[i]),input$selectTaxonomy]))
    }
    taxon_group
  }
})

family_group12_2 <- reactive({
  if (input$OmicTable == "OTUs"){
    taxon_group <- c()
    for (i in 1:nrow(corr_expr())){
      taxon_group <- c(taxon_group, (taxTable2()[which(taxTable2()[["rn"]] == rownames(corr_expr())[i]),input$selectTaxonomy]))
    }
    taxon_group
  }

})

family_group23_2 <- reactive({
  if (input$OmicTable == "OTUs"){
    taxon_group <- c()
    for (i in 1:ncol(corr_expr())){
      taxon_group <- c(taxon_group, (taxTable2()[which(taxTable2()[["rn"]] == colnames(corr_expr23())[i]),input$selectTaxonomy]))
    }
    taxon_group
  }
})

family_group23_3 <- reactive({
  if (input$OmicTable3 == "OTUs"){
    taxon_group <- c()
    for (i in 1:nrow(corr_expr())){
      taxon_group <- c(taxon_group, (taxTable2()[which(taxTable2()[["rn"]] == rownames(corr_expr23())[i]),input$selectTaxonomy]))
    }
    taxon_group
  }

})

family_group13_1 <- reactive({
  if (input$CountingT1 == "OTUs"){
    taxon_group <- c()
    for (i in 1:ncol(corr_expr())){
      taxon_group <- c(taxon_group, (taxTable2()[which(taxTable2()[["rn"]] == colnames(corr_expr13())[i]),input$selectTaxonomy]))
    }
    taxon_group
  }
})

family_group23_3 <- reactive({
  if (input$OmicTable3 == "OTUs"){
    taxon_group <- c()
    for (i in 1:nrow(corr_expr())){
      taxon_group <- c(taxon_group, (taxTable2()[which(taxTable2()[["rn"]] == rownames(corr_expr13())[i]),input$selectTaxonomy]))
    }
    taxon_group
  }

})

#### BI-PARTITE NETWORK
adjacency_expr12 <- reactive({
  corr_expr_1 <- corr_expr()
  #corr_expr_1[which(corr_expr_1 < 0.30 && corr_expr_1 > -0.30),] <- 0
  corr_expr_1[abs(corr_expr_1) < 0.5] <- 0
  #corr_expr_1[(corr_expr_1 < 0.70 && corr_expr_1 > 0.50) || (corr_expr_1 > -0.70 && corr_expr_1 < -0.50)] <- 1
  corr_expr_1[abs(corr_expr_1)>=0.5] <- 1
  #corr_expr_1 <- t(corr_expr_1)
  # if (nchar(rownames(corr_expr_1)[5]) >20){
  #   for (row in 1:nrow(corr_expr_1)){
  #
  #     rownames(corr_expr_1)[row] <- paste(input$CountingT1, "_D1", row, sep="")
  #   }
  # }
  # if (nchar(colnames(corr_expr_1)[5]) >20){
  #   for (col in 1:ncol(corr_expr_1)){
  #
  #     colnames(corr_expr_1)[col] <- paste(input$OmicTable, "_D2", row, sep="")
  #   }
  # }
  adjacency_expr <- corr_expr_1
  adjacency_expr
})


adjacency_expr13 <- reactive({
  corr_expr_1 <- corr_expr13()
  #corr_expr_1[which(corr_expr_1 < 0.30 && corr_expr_1 > -0.30),] <- 0
  corr_expr_1[abs(corr_expr_1) < 0.5] <- 0
  #corr_expr_1[(corr_expr_1 < 0.70 && corr_expr_1 > 0.50) || (corr_expr_1 > -0.70 && corr_expr_1 < -0.50)] <- 1
  corr_expr_1[abs(corr_expr_1)>=0.5] <- 1
  #corr_expr_1 <- t(corr_expr_1)
  # if (nchar(rownames(corr_expr_1)[5]) >20){
  #   for (row in 1:nrow(corr_expr_1)){
  #
  #     rownames(corr_expr_1)[row] <- paste(input$CountingT1, "_D1", row, sep="")
  #   }
  # }
  # if (nchar(colnames(corr_expr_1)[5]) >20){
  #   for (col in 1:ncol(corr_expr_1)){
  #
  #     colnames(corr_expr_1)[col] <- paste(input$OmicTable3, "_D3", row, sep="")
  #   }
  # }
  adjacency_expr <- corr_expr_1
  adjacency_expr
})

adjacency_expr23 <- reactive({
  corr_expr_1 <- corr_expr23()
  #corr_expr_1[which(corr_expr_1 < 0.30 && corr_expr_1 > -0.30),] <- 0
  corr_expr_1[abs(corr_expr_1) < 0.5] <- 0
  #corr_expr_1[(corr_expr_1 < 0.70 && corr_expr_1 > 0.50) || (corr_expr_1 > -0.70 && corr_expr_1 < -0.50)] <- 1
  corr_expr_1[abs(corr_expr_1)>=0.5] <- 1
  #corr_expr_1 <- t(corr_expr_1)
  # if (nchar(rownames(corr_expr_1)[5]) >20){
  #   for (row in 1:nrow(corr_expr_1)){
  #
  #     rownames(corr_expr_1)[row] <- paste(input$OmicTable, "_D2", row, sep="")
  #   }
  # }
  # if (nchar(colnames(corr_expr_1)[5]) >20){
  #   for (col in 1:ncol(corr_expr_1)){
  #
  #     colnames(corr_expr_1)[col] <- paste(input$OmicTable3, "_D3", row, sep="")
  #   }
  # }
  adjacency_expr <- corr_expr_1
  adjacency_expr
})

bip12 <- reactive({

  bip = network::network(adjacency_expr12(),
                         matrix.type = "bipartite",
                         ignore.eval = FALSE,
                         names.eval = "weights")


  bip
})

bip23 <- reactive({

  bip = network::network(adjacency_expr23(),
                         matrix.type = "bipartite",
                         ignore.eval = FALSE,
                         names.eval = "weights")


  bip
})

bip13 <- reactive({

  bip = network::network(adjacency_expr13(),
                         matrix.type = "bipartite",
                         ignore.eval = FALSE,
                         names.eval = "weights")


  bip
})


#### INTERFACE VARIABLES ####


# output$SelectModule1 <- renderUI({
#   selectInput("selectModule1_p5",
#               label= "Choose a module color: ",
#               choices = c("General", substr(names(selectedMEs()), 3, 40)),
#               selected = "Individual Variables")
# })
#
# output$SelectModule3 <- renderUI({
#   selectInput("selectModule3_p5",
#               label= "Choose a module color: ",
#               choices = c("General", substr(names(selectedMEs3()), 3, 40)),
#               selected = "Individual Variables")
# })

output$SelectModule12 <- renderUI({
  selectInput("selectModule12_p5",
              label= "Choose a module color: ",
              choices = c("Individual_Variables", substr(names(selectedMEs()), 3, 40)),
              selected = NULL)
})

output$SelectModule13 <- renderUI({
  selectInput("selectModule13_p5",
              label= "Choose a module color: ",
              choices = c("Individual_Variables", substr(names(selectedMEs()), 3, 40)),
              selected = NULL)
})

# output$SelectModule2 <- renderUI({
#   selectInput("selectModule2_p5",
#               label= "Choose a module color: ",
#               choices = c("General", substr(names(selectedMEs2()), 3, 40)),
#               selected = "General")
# })

output$SelectModule21 <- renderUI({
  selectInput("selectModule21_p5",
              label= "Choose a module color: ",
              choices = c("Individual_Variables", substr(names(selectedMEs2()), 3, 40)),
              selected = NULL)
})

output$SelectModule23 <- renderUI({
  selectInput("selectModule23_p5",
              label= "Choose a module color: ",
              choices = c("Individual_Variables", substr(names(selectedMEs2()), 3, 40)),
              selected = NULL)
})

output$SelectModule31 <- renderUI({
  selectInput("selectModule31_p5",
              label= "Choose a module color: ",
              choices = c("Individual_Variables", substr(names(selectedMEs3()), 3, 40)),
              selected = NULL)
})

output$SelectModule32 <- renderUI({
  selectInput("selectModule32_p5",
              label= "Choose a module color: ",
              choices = c("Individual_Variables", substr(names(selectedMEs3()), 3, 40)),
              selected = NULL)
})



output$SelectVariable2 <- renderUI({
  selectInput("SelectVariable2",
              label = "Choose a variable: ",
              choices = colnames(sampleAnnot_2()),
              selected = colnames(sampleAnnot_2())[2])
})

output$traitHive <- renderUI({
  selectInput("traitHive",
              label = "Choose a variable: ",
              choices = c(colnames(sampleAnnot_2())),
              selected = colnames(sampleAnnot_2())[2])
})

output$selectVariables <- renderUI({
  selectInput("selectVariables",
              label = "Choose variables to plot",
              choices = colnames(exprDatSec_2()),
              multiple = TRUE)
})

output$selectVariablesExprDat <- renderUI({
  selectInput("selectVariablesExprDat",
              label = "Choose variables to plot",
              choices = colnames(exprDat_3()),
              multiple = TRUE)
})

output$selectVariables31 <- renderUI({
  selectInput("selectVariables31",
              label = "Choose variables to plot",
              choices = colnames(exprDatTer_3()),
              multiple = TRUE)
})

output$selectVariablesExprDat13 <- renderUI({
  selectInput("selectVariablesExprDat13",
              label = "Choose variables to plot",
              choices = colnames(exprDat_3()),
              multiple = TRUE)
})

output$selectVariables32 <- renderUI({
  selectInput("selectVariables32",
              label = "Choose variables to plot",
              choices = colnames(exprDatTer_2()),
              multiple = TRUE)
})

output$selectVariablesExprDat23 <- renderUI({
  selectInput("selectVariablesExprDat23",
              label = "Choose variables to plot",
              choices = colnames(exprDatSec_2()),
              multiple = TRUE)
})

output$SelectTaxonomy <- renderUI({
  taxonomy <- colnames(taxTable1())
  taxonomy <- taxonomy[2:ncol(taxTable1())]
  selectInput("selectTaxonomy",
              label = "Select a taxnomic level",
              choices = taxonomy,
              selected = taxonomy[5])
})

#### PLOT OUTPUTS ####

output$coinertia <- renderPlot({
  if (is.null(selected_sampleInfo_all()))
    return(NULL)
  if (!input$Omic3){
    if (!is.numeric(sampleAnnot_2()[,input$SelectVariable2])){
      if (input$ShowDrivers){
        coinertia_plot <-
          ggplot(data = selected_sampleInfo_all()) +
          geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[7]), size = 3, shape = 15) +
          geom_point(aes_string(x = colnames(selected_sampleInfo_all())[3], y = colnames(selected_sampleInfo_all())[4], col = colnames(selected_sampleInfo_all())[7]), shape = 16, size = 3) +
          scale_shape(solid = TRUE) +
          theme_bw() +
          geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[7], xend = colnames(selected_sampleInfo_all())[3], yend = colnames(selected_sampleInfo_all())[4])) +
          geom_label_repel(data = selected_axis1_drivers(),
                           aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
                           size = 4, segment.size = 0.5,
                           label.padding = unit(0.1, "lines"), label.size = 0.5) +
          labs(x = sprintf("Axis1 [%s%% Variance]",
                           100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
               y = sprintf("Axis2 [%s%% Variance]",
                           100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)))
      }else{
        coinertia_plot <-
          ggplot(data = selected_sampleInfo_all()) +
          geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[7]), size = 3, shape = 15) +
          geom_point(aes_string(x = colnames(selected_sampleInfo_all())[3], y = colnames(selected_sampleInfo_all())[4], col = colnames(selected_sampleInfo_all())[7]), shape = 16, size = 3) +
          scale_shape(solid = TRUE) +
          theme_bw() +
          geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[7], xend = colnames(selected_sampleInfo_all())[3], yend = colnames(selected_sampleInfo_all())[4])) +
          # geom_label_repel(data = selected_axis1_drivers(),
          #                  aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
          #                  size = 4, segment.size = 0.5,
          #                  label.padding = unit(0.1, "lines"), label.size = 0.5) +
          labs(x = sprintf("Axis1 [%s%% Variance]",
                           100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
               y = sprintf("Axis2 [%s%% Variance]",
                           100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)))
      }
    }else{
      if (input$ShowDrivers){
        coinertia_plot <-
          ggplot(data = selected_sampleInfo_all()) +
          geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[7]), size = 3, shape = 15) +
          geom_point(aes_string(x = colnames(selected_sampleInfo_all())[3], y = colnames(selected_sampleInfo_all())[4], col = colnames(selected_sampleInfo_all())[7]), shape = 16, size = 3) +
          scale_shape(solid = TRUE) +
          theme_bw() +
          scale_colour_gradientn(colours = terrain.colors(10)) +
          geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[7], xend = colnames(selected_sampleInfo_all())[3], yend = colnames(selected_sampleInfo_all())[4])) +
          geom_label_repel(data = selected_axis1_drivers(),
                           aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
                           size = 4, segment.size = 0.5,
                           label.padding = unit(0.1, "lines"), label.size = 0.5) +
          labs(x = sprintf("Axis1 [%s%% Variance]",
                           100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
               y = sprintf("Axis2 [%s%% Variance]",
                           100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)))
      }else{
        coinertia_plot <-
          ggplot(data = selected_sampleInfo_all()) +
          geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[7]), size = 3, shape = 15) +
          geom_point(aes_string(x = colnames(selected_sampleInfo_all())[3], y = colnames(selected_sampleInfo_all())[4], col = colnames(selected_sampleInfo_all())[7]), shape = 16, size = 3) +
          scale_shape(solid = TRUE) +
          theme_bw() +
          scale_colour_gradientn(colours = terrain.colors(10)) +
          geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[7], xend = colnames(selected_sampleInfo_all())[3], yend = colnames(selected_sampleInfo_all())[4])) +
          # geom_label_repel(data = selected_axis1_drivers(),
          #                  aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
          #                  size = 4, segment.size = 0.5,
          #                  label.padding = unit(0.1, "lines"), label.size = 0.5) +
          labs(x = sprintf("Axis1 [%s%% Variance]",
                           100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
               y = sprintf("Axis2 [%s%% Variance]",
                           100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)))
      }

    }

    bivariate_plot <-
      ggplot(data = bivariate_plot(),
             aes_string(x = colnames(bivariate_plot())[1],
                        y = colnames(bivariate_plot())[2],
                        label = "sampleName")) +
      theme_bw() +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, col = "red") +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      geom_text(vjust = 0, nudge_y = 0.02) +
      annotate("text", label = "maximal correlation", x = -0.5, y =-0.5, color = "red", size = 5)
    print(plot_grid(coinertia_plot,LegendDF(),bivariate_plot, ncol = 2, rel_widths = c(1, .3)))
  }else{


    all_sample_info <- selected_sampleInfo_all()
    mco <- mcoia()
    axis1_drivers <- selected_axis1_drivers()
    if (!is.numeric(sampleAnnot_2()[,input$SelectVariable2])){
      if (input$ShowDrivers){
        coinertia_plot <- ggplot(data = all_sample_info) +
          geom_point(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], col = colnames(all_sample_info)[9]), size = 3, shape = 15) +
          geom_point(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], col = colnames(all_sample_info)[9]), shape = 16, size = 3) +
          geom_point(aes_string(x = colnames(all_sample_info)[5], y = colnames(all_sample_info)[6], col = colnames(all_sample_info)[9]), shape = 17, size = 3) +
          geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[3], yend = colnames(all_sample_info)[4])) +
          geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
          geom_segment(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
          geom_label_repel(data = axis1_drivers,
                           aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
                           size = 4, segment.size = 0.5,
                           label.padding = unit(0.1, "lines"), label.size = 0.5) +
          theme_bw() +
          labs(x = sprintf("Axis1 [%s%% Variance]",
                           100 * round(mco[["mcoa"]]$pseudoeig[1] / sum(mco[["mcoa"]]$pseudoeig), 2)),
               y = sprintf("Axis2 [%s%% Variance]",
                           100 * round(mco[["mcoa"]]$pseudoeig[2] / sum(mco[["mcoa"]]$pseudoeig), 2)))
      }else{
        coinertia_plot <- ggplot(data = all_sample_info) +
          geom_point(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], col = colnames(all_sample_info)[9]), size = 3, shape = 15) +
          geom_point(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], col = colnames(all_sample_info)[9]), shape = 16, size = 3) +
          geom_point(aes_string(x = colnames(all_sample_info)[5], y = colnames(all_sample_info)[6], col = colnames(all_sample_info)[9]), shape = 17, size = 3) +
          geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[3], yend = colnames(all_sample_info)[4])) +
          geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
          geom_segment(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
          # geom_label_repel(data = axis1_drivers,
          #                  aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
          #                  size = 4, segment.size = 0.5,
          #                  label.padding = unit(0.1, "lines"), label.size = 0.5) +
          theme_bw() +
          labs(x = sprintf("Axis1 [%s%% Variance]",
                           100 * round(mco[["mcoa"]]$pseudoeig[1] / sum(mco[["mcoa"]]$pseudoeig), 2)),
               y = sprintf("Axis2 [%s%% Variance]",
                           100 * round(mco[["mcoa"]]$pseudoeig[2] / sum(mco[["mcoa"]]$pseudoeig), 2)))
      }

    }else{
      if (input$ShowDrivers){
        coinertia_plot <- ggplot(data = all_sample_info) +
          geom_point(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], col = colnames(all_sample_info)[9]), size = 3, shape = 15) +
          geom_point(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], col = colnames(all_sample_info)[9]), shape = 16, size = 3) +
          geom_point(aes_string(x = colnames(all_sample_info)[5], y = colnames(all_sample_info)[6], col = colnames(all_sample_info)[9]), shape = 17, size = 3) +
          geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[3], yend = colnames(all_sample_info)[4])) +
          geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
          geom_segment(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
          geom_label_repel(data = axis1_drivers,
                           aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
                           size = 4, segment.size = 0.5,
                           label.padding = unit(0.1, "lines"), label.size = 0.5) +
          theme_bw() +
          scale_colour_gradientn(colours = terrain.colors(10)) +
          labs(x = sprintf("Axis1 [%s%% Variance]",
                           100 * round(mco[["mcoa"]]$pseudoeig[1] / sum(mco[["mcoa"]]$pseudoeig), 2)),
               y = sprintf("Axis2 [%s%% Variance]",
                           100 * round(mco[["mcoa"]]$pseudoeig[2] / sum(mco[["mcoa"]]$pseudoeig), 2)))
      }else{

        coinertia_plot <- ggplot(data = all_sample_info) +
          geom_point(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], col = colnames(all_sample_info)[9]), size = 3, shape = 15) +
          geom_point(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], col = colnames(all_sample_info)[9]), shape = 16, size = 3) +
          geom_point(aes_string(x = colnames(all_sample_info)[5], y = colnames(all_sample_info)[6], col = colnames(all_sample_info)[9]), shape = 17, size = 3) +
          geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[3], yend = colnames(all_sample_info)[4])) +
          geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
          geom_segment(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
          # geom_label_repel(data = axis1_drivers,
          #                  aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
          #                  size = 4, segment.size = 0.5,
          #                  label.padding = unit(0.1, "lines"), label.size = 0.5) +
          theme_bw() +
          scale_colour_gradientn(colours = terrain.colors(10)) +
          labs(x = sprintf("Axis1 [%s%% Variance]",
                           100 * round(mco[["mcoa"]]$pseudoeig[1] / sum(mco[["mcoa"]]$pseudoeig), 2)),
               y = sprintf("Axis2 [%s%% Variance]",
                           100 * round(mco[["mcoa"]]$pseudoeig[2] / sum(mco[["mcoa"]]$pseudoeig), 2)))
      }

    }

    print(plot_grid(coinertia_plot,LegendDF(), ncol = 2, rel_widths = c(1, .3)))
  }



})

output$coinertia_bivariate <- renderPlot({
  if (is.null(distTree_coia()))
    return(NULL)
  validate(
    need(!input$Omic3, 'Only for multi-omics analysis with two datasets')
  )
  dendrogramme_coia <-
    ggdendrogram(distTree_coia())
  plot_dist_coia <-
    ggplot(data = distance_mcoia_df_supp_info(), aes(x= sampleName, y = distance_coia.v, col = color)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot_hist_coia <-
    ggplot(distance_mcoia_df_supp_info(), aes(distance_coia.v, fill = color)) + geom_histogram(binwidth = 0.1)
  print(plot_grid(dendrogramme_coia,
                  plot_dist_coia,
                  plot_hist_coia,
                  labels = c("A", "B", "C"),
                  nrow = 2,
                  align = "h"))

})

output$RV <- renderPrint({
  if (is.null(RV_score()))
    return(NULL)
  RV_score()
})

output$RV2 <- renderPrint({
  if (is.null(RV_score()))
    return(NULL)
  RV_score()
})

output$PA <- renderPlot({
  if (is.null(sample_info_PA()))
    return(NULL)
  validate(
    need(!input$Omic3, 'Only for multi-omics analysis with two datasets')
  )
  ggplot(data = sample_info_PA()) + geom_point(aes_string(x = 'Axis1', y= 'Axis2', col = input$SelectVariable2), size = 3) +
    geom_point(aes_string(x = 'Axis.1', y= 'Axis.2', col = input$SelectVariable2), size = 3, shape = 5) +
    geom_segment(aes_string(x = 'Axis1', y = 'Axis2', xend = 'Axis.1', yend = 'Axis.2', colour = input$SelectVariable2))
})


output$HIVE_MEs <- renderPlot({
  # if (input$traitHive == "Keystone"){
  #   myHive <- hive_data_keystone()
  # }else{
    myHive <- hive_data()
  # }
  if (is.null(myHive))
    return(NULL)

  myHive$nodes$lab <- as.character(myHive$nodes$lab)
  myHive$nodes$color <- as.character(myHive$nodes$color)
  myHive$edges$color <- as.character(myHive$edges$color)
  if (input$Omic3){
    plotMyHive(myHive, bkgnd = "white", axLabs = c(paste(input$CountingT1, "D1", sep = "_"), paste( input$OmicTable, "D2", sep = "_"), paste(input$OmicTable3, "D3", sep = "_")), axLab.gpar = gpar(col = "#636363"))
  }else{
    plotMyHive(myHive, ch = 2, bkgnd = "white", axLabs = c(paste(input$CountingT1, "D1", sep = "_"), paste( input$OmicTable, "D2", sep = "_")), axLab.gpar = gpar(col = "#636363"), rot = c(-90, 90))

  }

})

# output$NETWORK_MEs <- renderPlot({
#   myHive <- hive_data_keystone()
#   Keystone_net <- keystoneMyNodes(nodes=myHive$nodes, edges= myHive$edges, returnNetwork= TRUE)
#   plot(Keystone_net, vertex.label=V(Keystone_net)$lab, vertex.label.cex = .7)
# })

output$HEATMAP_MEs12 <- renderIheatmap({
  if (is.null(corr_MEs_D1D2()))
    return(NULL)
  main_heatmap(corr_MEs_D1D2(), name = "Correlation between MEs") %>%
    add_row_labels() %>%
    add_col_labels() %>%
    add_col_clustering() %>%
    add_row_clustering() %>%
    add_col_title("Modules of the first dataset") %>%
    add_row_title("Modules of the second dataset")

})

output$HEATMAP_MEs23 <- renderIheatmap({
  if (is.null(corr_MEs_D2D3()))
    return(NULL)
  main_heatmap(corr_MEs_D2D3(), name = "Correlation between MEs") %>%
    add_row_labels() %>%
    add_col_labels() %>%
    add_col_clustering() %>%
    add_row_clustering() %>%
    add_col_title("Modules of the first dataset") %>%
    add_row_title("Modules of the second dataset")

})

output$HEATMAP_MEs13 <- renderIheatmap({
  if (is.null(corr_MEs_D1D3()))
    return(NULL)
  main_heatmap(corr_MEs_D1D3(), name = "Correlation between MEs") %>%
    add_row_labels() %>%
    add_col_labels() %>%
    add_col_clustering() %>%
    add_row_clustering() %>%
    add_col_title("Modules of the first dataset") %>%
    add_row_title("Modules of the second dataset")

})

output$HEATMAP <- renderIheatmap({
  if (is.null(corr_expr()))
    return(NULL)
  if (input$addLabels12){
    hm <- main_heatmap(corr_expr(), name = "Correlation") %>%
      add_row_labels() %>%
      add_col_labels() %>%
      add_col_clustering() %>%
      add_row_clustering() %>%
      add_col_title("Variables of the first datasets in each module") %>%
      add_row_title("Variables of the second dataset")
  }else{
    hm <- main_heatmap(corr_expr(), name = "Correlation") %>%
      add_col_clustering() %>%
      add_row_clustering() %>%
      add_col_title("Variables of the first datasets in each module") %>%
      add_row_title("Variables of the second dataset")
  }
  hm

})

output$HEATMAP23 <- renderIheatmap({
  if (is.null(corr_expr23()))
    return(NULL)
  if (input$addLabels23){
    hm <- main_heatmap(corr_expr23(), name = "Correlation") %>%
      add_row_labels() %>%
      add_col_labels() %>%
      add_col_clustering() %>%
      add_row_clustering() %>%
      add_col_title("Variables of the first datasets in each module") %>%
      add_row_title("Variables of the second dataset")
  }else{
    hm <- main_heatmap(corr_expr23(), name = "Correlation") %>%
      add_col_clustering() %>%
      add_row_clustering() %>%
      add_col_title("Variables of the first datasets in each module") %>%
      add_row_title("Variables of the second dataset")
  }
  hm
})

output$HEATMAP13 <- renderIheatmap({
  if (is.null(corr_expr13()))
    return(NULL)
  if (input$addLabels13){
    hm <- main_heatmap(corr_expr13(), name = "Correlation") %>%
      add_row_labels() %>%
      add_col_labels() %>%
      add_col_clustering() %>%
      add_row_clustering() %>%
      add_col_title("Variables of the first datasets in each module") %>%
      add_row_title("Variables of the second dataset")
  }else{
    hm <- main_heatmap(corr_expr13(), name = "Correlation") %>%
      add_col_clustering() %>%
      add_row_clustering() %>%
      add_col_title("Variables of the first datasets in each module") %>%
      add_row_title("Variables of the second dataset")
  }
})

output$bipartite_network12 <- renderPlot({
  if (is.null(bip12()))
    return(NULL)
  col = c("actor" = "red", "event" = "blue")
  ggnet2(bip12(), color = "mode", palette = col, label = TRUE, layout.exp = 0.25, shape= "mode") +
    theme(legend.position="none")

})

output$bipartite_network23 <- renderPlot({
  if (is.null(bip23()))
    return(NULL)
  col = c("actor" = "red", "event" = "blue")
  ggnet2(bip23(), color = "mode", palette = col, label = TRUE, layout.exp = 0.25, shape= "mode") +
    theme(legend.position="none")

})

output$bipartite_network13 <- renderPlot({
  if (is.null(bip13()))
    return(NULL)
  col = c("actor" = "red", "event" = "blue")
  ggnet2(bip13(), color = "mode", palette = col, label = TRUE, layout.exp = 0.25, shape= "mode") +
    theme(legend.position="none")

})

#### DOWNLOADS ####

output$Download_Coinertia_drivers <- downloadHandler(
  filename = function(){
    if (input$Omic3){
      paste(input$selectModule1_p5, "_", input$selectModule2_p5, "_", input$selectModule3_p5, "coinertia_drivers.csv", sep="")
    }else{
      paste(input$selectModule1_p5, "_", input$selectModule2_p5, "coinertia_drivers.csv", sep="")
    }
  },
  content = function (filename){
    write.csv(as.data.frame(mcoia()[["mcoa"]]$Tco), filename)
  }
)


output$Download_Multivariate_Analysis2 <- downloadHandler(
  filename = function(){
    paste("HIVEPLOT.zip")
  },
  content = function (filename){
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    myHive <- hive_data()
    validate(
      need(nrow(myHive$edges) != 0, "No significant interaction between modules and external trait")
    )
    myHive$nodes$lab <- as.character(myHive$nodes$lab)
    myHive$nodes$color <- as.character(myHive$nodes$color)
    myHive$edges$color <- as.character(myHive$edges$color)
    if (input$Omic3){
      if (input$pdf_or_svg_p5_2 == "pdf"){
        pdf(paste("hive.pdf", sep = ""), width = input$widthPDF2, height = input$heightPDF2)
        plotMyHive(myHive, bkgnd = "white", axLabs = c(paste(input$CountingT1, "D1", sep = "_"), paste( input$OmicTable, "D2", sep = "_"), paste(input$OmicTable3, "D3", sep = "_")), axLab.gpar = gpar(col = "#636363"))
        dev.off()
        fs <- c(fs, paste("hive.pdf"))
      }else{
        svg(paste("hive.svg", sep = ""), width = input$widthPDF2, height = input$heightPDF2)
        plotMyHive(myHive, bkgnd = "white", axLabs = c(paste(input$CountingT1, "D1", sep = "_"), paste( input$OmicTable, "D2", sep = "_"), paste(input$OmicTable3, "D3", sep = "_")), axLab.gpar = gpar(col = "#636363"))
        dev.off()
        fs <- c(fs, paste("hive.svg"))
      }
    }else{
      if (input$pdf_or_svg_p5_2 == "pdf"){
        pdf(paste("hive.pdf", sep = ""), width = input$widthPDF2, height = input$heightPDF2)
        plotMyHive(myHive, ch = 2, bkgnd = "white", axLabs = c(paste(input$CountingT1, "D1", sep = "_"), paste( input$OmicTable, "D2", sep = "_")), axLab.gpar = gpar(col = "#636363"), rot = c(-90, 90))
        dev.off()
        fs <- c(fs, paste("hive.pdf"))
      }else{
        svg(paste("hive.svg", sep = ""), width = input$widthPDF2, height = input$heightPDF2)
        plotMyHive(myHive, ch = 2, bkgnd = "white", axLabs = c(paste(input$CountingT1, "D1", sep = "_"), paste( input$OmicTable, "D2", sep = "_")), axLab.gpar = gpar(col = "#636363"), rot = c(-90, 90))
        dev.off()
        fs <- c(fs, paste("hive.svg"))
      }
    }
    write.csv(myHive$edges, "Hive_edges.csv")
    fs <- c(fs, "Hive_edges.csv")
    write.csv(myHive$nodes, "Hive_nodes.csv")
    fs <- c(fs, "Hive_nodes.csv")

    zip(zipfile=filename, files=fs)


  },
  contentType = "application/zip")

output$Download_Multivariate_Analysis3 <- downloadHandler(
  filename = function(){
    paste("HEATMAP12.zip")
  },
  content = function (filename){
    webshot::install_phantomjs()
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    hm <- main_heatmap(corr_MEs_D1D2(), name = "Correlation between MEs") %>%
      add_row_labels() %>%
      add_col_labels() %>%
      add_col_clustering() %>%
      add_row_clustering() %>%
      add_col_title("Modules of the first dataset") %>%
      add_row_title("Modules of the second dataset")

    if (input$pdf_or_svg_p5_3 == "pdf"){
      save_iheatmap(hm, "heatmap12.pdf" )
      fs <- c(fs, paste("heatmap12.pdf"))
    }else{
      save_iheatmap(hm, "heatmap12.png" )
      fs <- c(fs, paste("heatmap12.png"))
    }
    write.csv(corr_MEs_D1D2(), "Corr_D1D2.csv")
    fs <- c(fs, "Corr_D1D2.csv")
    write.csv(p_val_MEs_D1D2(), "PValues_D1D2.csv")
    fs <- c(fs, "PValues_D1D2.csv")

    zip(zipfile=filename, files=fs)


  },
  contentType = "application/zip")


output$Download_Multivariate_Analysis4 <- downloadHandler(
  filename = function(){
    paste("HEATMAP23.zip")
  },
  content = function (filename){
    webshot::install_phantomjs()
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    hm <- main_heatmap(corr_MEs_D2D3(), name = "Correlation between MEs") %>%
      add_row_labels() %>%
      add_col_labels() %>%
      add_col_clustering() %>%
      add_row_clustering() %>%
      add_col_title("Modules of the second dataset") %>%
      add_row_title("Modules of the third dataset")

    if (input$pdf_or_svg_p5_4 == "pdf"){
      save_iheatmap(hm, "heatmap23.pdf" )
      fs <- c(fs, paste("heatmap23.pdf"))
    }else{
      save_iheatmap(hm, "heatmap23.png" )
      fs <- c(fs, paste("heatmap23.png"))
    }
    write.csv(corr_MEs_D2D3(), "Corr_D2D3.csv")
    fs <- c(fs, "Corr_D2D3.csv")
    write.csv(p_val_MEs_D2D3(), "PValues_D2D3.csv")
    fs <- c(fs, "PValues_D2D3.csv")

    zip(zipfile=filename, files=fs)


  },
  contentType = "application/zip")


output$Download_Multivariate_Analysis5 <- downloadHandler(
  filename = function(){
    paste("HEATMAP13.zip")
  },
  content = function (filename){
    webshot::install_phantomjs()
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    hm <- main_heatmap(corr_MEs_D1D3(), name = "Correlation between MEs") %>%
      add_row_labels() %>%
      add_col_labels() %>%
      add_col_clustering() %>%
      add_row_clustering() %>%
      add_col_title("Modules of the first dataset") %>%
      add_row_title("Modules of the third dataset")

    if (input$pdf_or_svg_p5_5 == "pdf"){
      save_iheatmap(hm, "heatmap13.pdf" )
      fs <- c(fs, paste("heatmap13.pdf"))
    }else{
      save_iheatmap(hm, "heatmap13.png" )
      fs <- c(fs, paste("heatmap13.png"))
    }
    write.csv(corr_MEs_D1D3(), "Corr_D1D3.csv")
    fs <- c(fs, "Corr_D1D3.csv")
    write.csv(p_val_MEs_D1D3(), "PValues_D1D3.csv")
    fs <- c(fs, "PValues_D1D3.csv")

    zip(zipfile=filename, files=fs)


  },
  contentType = "application/zip")


output$Download_Multivariate_Analysis6 <- downloadHandler(
  filename = function(){
    paste("HEATMAP_AND_BIPARTITE.zip")
  },
  content = function (filename){
    webshot::install_phantomjs()
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    if (input$addLabels12){
      hm <- main_heatmap(corr_expr(), name = "Correlation") %>%
        add_row_labels() %>%
        add_col_labels() %>%
        add_col_clustering() %>%
        add_row_clustering() %>%
        add_col_title("Variables of the first datasets in each module") %>%
        add_row_title("Variables of the second dataset")
    }else{
      hm <- main_heatmap(corr_expr(), name = "Correlation") %>%
        add_col_clustering() %>%
        add_row_clustering() %>%
        add_col_title("Variables of the first datasets in each module") %>%
        add_row_title("Variables of the second dataset")
    }


    if (input$pdf_or_svg_p5_6 == "pdf"){
      save_iheatmap(hm, paste("heatmap12_",input$selectModule12_p5,"_", input$selectModule21_p5, ".pdf", sep = "") )
      fs <- c(fs, paste("heatmap12_",input$selectModule12_p5,"_", input$selectModule21_p5, ".pdf", sep = ""))
      pdf(paste("bipartite_",input$selectModule12_p5,"_", input$selectModule21_p5, ".pdf", sep = ""), height = input$heightPDF6, width = input$widthPDF6)
      col = c("actor" = "red", "event" = "blue")
      print(ggnet2(bip12(), color = "mode", palette = col, label = TRUE, layout.exp = 0.25, shape= "mode") +
              theme(legend.position="none")    )
      dev.off()
      fs <- c(fs, paste("bipartite_",input$selectModule12_p5,"_", input$selectModule21_p5, ".pdf", sep = ""))
    }else{
      save_iheatmap(hm, paste("heatmap12_",input$selectModule12_p5,"_", input$selectModule21_p5, ".png", sep = "") )
      fs <- c(fs, paste("heatmap12_",input$selectModule12_p5,"_", input$selectModule21_p5, ".png", sep = ""))
      svg(paste("bipartite_",input$selectModule12_p5,"_", input$selectModule21_p5, ".svg", sep = ""), height = input$heightPDF6, width = input$widthPDF6)
      col = c("actor" = "red", "event" = "blue")
      print(ggnet2(bip12(), color = "mode", palette = col, label = TRUE, layout.exp = 0.25, shape= "mode") +
              theme(legend.position="none")    )
      dev.off()
      fs <- c(fs, paste("bipartite_",input$selectModule12_p5,"_", input$selectModule21_p5, ".svg", sep = ""))
    }
    write.csv(corr_expr(), "Corr_expr_D1D2.csv")
    fs <- c(fs, "Corr_expr_D1D2.csv")


    zip(zipfile=filename, files=fs)


  },
  contentType = "application/zip")

output$Download_Multivariate_Analysis7 <- downloadHandler(
  filename = function(){
    paste("HEATMAP_AND_BIPARTITE.zip")
  },
  content = function (filename){
    webshot::install_phantomjs()
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    if (input$addLabels23){
      hm <- main_heatmap(corr_expr23(), name = "Correlation") %>%
        add_row_labels() %>%
        add_col_labels() %>%
        add_col_clustering() %>%
        add_row_clustering() %>%
        add_col_title("Variables of the first datasets in each module") %>%
        add_row_title("Variables of the second dataset")
    }else{
      hm <- main_heatmap(corr_expr23(), name = "Correlation") %>%
        add_col_clustering() %>%
        add_row_clustering() %>%
        add_col_title("Variables of the first datasets in each module") %>%
        add_row_title("Variables of the second dataset")
    }

    if (input$pdf_or_svg_p5_7 == "pdf"){
      save_iheatmap(hm, paste("heatmap23_",input$selectModule23_p5,"_", input$selectModule32_p5, ".pdf", sep = "") )
      fs <- c(fs, paste("heatmap23_",input$selectModule23_p5,"_", input$selectModule32_p5, ".pdf", sep = ""))
      pdf(paste("bipartite_",input$selectModule23_p5,"_", input$selectModule32_p5, ".pdf", sep = ""), height = input$heightPDF7, width = input$widthPDF7)
      col = c("actor" = "red", "event" = "blue")
      print(ggnet2(bip23(), color = "mode", palette = col, label = TRUE, layout.exp = 0.25, shape= "mode") +
              theme(legend.position="none")    )
      dev.off()
      fs <- c(fs, paste("bipartite_",input$selectModule23_p5,"_", input$selectModule32_p5, ".pdf", sep = ""))
    }else{
      save_iheatmap(hm, paste("heatmap23_",input$selectModule23_p5,"_", input$selectModule32_p5, ".png", sep = "") )
      fs <- c(fs, paste("heatmap23_",input$selectModule23_p5,"_", input$selectModule32_p5, ".png", sep = ""))
      svg(paste("bipartite_",input$selectModule23_p5,"_", input$selectModule32_p5, ".svg", sep = ""), height = input$heightPDF7, width = input$widthPDF7)
      col = c("actor" = "red", "event" = "blue")
      print(ggnet2(bip23(), color = "mode", palette = col, label = TRUE, layout.exp = 0.25, shape= "mode") +
              theme(legend.position="none")    )
      dev.off()
      fs <- c(fs, paste("bipartite_",input$selectModule23_p5,"_", input$selectModule32_p5, ".svg", sep = ""))
    }
    write.csv(corr_expr23(), "Corr_expr_D2D3.csv")
    fs <- c(fs, "Corr_expr_D2D3.csv")

    zip(zipfile=filename, files=fs)


  },
  contentType = "application/zip")


output$Download_Multivariate_Analysis8 <- downloadHandler(
  filename = function(){
    paste("HEATMAP_AND_BIPARTITE.zip")
  },
  content = function (filename){
    webshot::install_phantomjs()
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    if (input$addLabels13){
      hm <- main_heatmap(corr_expr13(), name = "Correlation") %>%
        add_row_labels() %>%
        add_col_labels() %>%
        add_col_clustering() %>%
        add_row_clustering() %>%
        add_col_title("Variables of the first datasets in each module") %>%
        add_row_title("Variables of the second dataset")
    }else{
      hm <- main_heatmap(corr_expr13(), name = "Correlation") %>%
        add_col_clustering() %>%
        add_row_clustering() %>%
        add_col_title("Variables of the first datasets in each module") %>%
        add_row_title("Variables of the second dataset")
    }

    if (input$pdf_or_svg_p5_8 == "pdf"){
      save_iheatmap(hm, paste("heatmap13_",input$selectModule13_p5,"_", input$selectModule31_p5, ".pdf", sep = "") )
      fs <- c(fs, paste("heatmap13_",input$selectModule13_p5,"_", input$selectModule31_p5, ".pdf", sep = ""))
      pdf(paste("bipartite_",input$selectModule13_p5,"_", input$selectModule31_p5, ".pdf", sep = ""), height = input$heightPDF8, width = input$widthPDF8)
      col = c("actor" = "red", "event" = "blue")
      print(ggnet2(bip13(), color = "mode", palette = col, label = TRUE, layout.exp = 0.25, shape= "mode") +
              theme(legend.position="none")    )
      dev.off()
      fs <- c(fs, paste("bipartite_",input$selectModule13_p5,"_", input$selectModule31_p5, ".pdf", sep = ""))
    }else{
      save_iheatmap(hm, paste("heatmap13_",input$selectModule13_p5,"_", input$selectModule31_p5, ".png", sep = "") )
      fs <- c(fs, paste("heatmap13_",input$selectModule13_p5,"_", input$selectModule31_p5, ".png", sep = ""))
      svg(paste("bipartite_",input$selectModule13_p5,"_", input$selectModule31_p5, ".svg", sep = ""), height = input$heightPDF8, width = input$widthPDF8)
      col = c("actor" = "red", "event" = "blue")
      print(ggnet2(bip13(), color = "mode", palette = col, label = TRUE, layout.exp = 0.25, shape= "mode") +
              theme(legend.position="none")  )
      dev.off()
      fs <- c(fs, paste("bipartite_",input$selectModule13_p5,"_", input$selectModule31_p5, ".svg", sep = ""))
    }
    write.csv(corr_expr(), "Corr_expr_D1D3.csv")
    fs <- c(fs, "Corr_expr_D1D3.csv")

    zip(zipfile=filename, files=fs)


  },
  contentType = "application/zip")

output$Download_Multivariate_Analysis <- downloadHandler(
  filename = function(){
    paste("ORDINATION.zip")
  },
  content = function (filename){
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    if (input$Omic3){
      if (input$pdf_or_svg_p5 == "pdf"){
        # pdf(paste("HIVEPLOT_", input$traitHive,".pdf", sep = ""), width = 10, height = 10)
        #
        # myHive <- hive_data()
        # myHive$nodes$lab <- as.character(myHive$nodes$lab)
        # myHive$nodes$color <- as.character(myHive$nodes$color)
        # myHive$edges$color <- as.character(myHive$edges$color)
        # plotMyHive(myHive, bkgnd = "white", axLabs = c(paste(input$CountingT1, "D1", sep = "_"), paste( input$OmicTable, "D2", sep = "_"), paste(input$OmicTable3, "D3", sep = "_")), axLab.gpar = gpar(col = "#636363"))
        # dev.off()
        #
        # fs <- c(fs,paste("HIVEPLOT_", input$traitHive,".pdf", sep = ""))


        pdf(paste("Co_inertia_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5, "_df3.",input$selectModule3_p5 , ".pdf", sep = ""), width = 10, height = 10)
        all_sample_info <- selected_sampleInfo_all()
        mco <- mcoia()
        axis1_drivers <- selected_axis1_drivers()
        if (!is.numeric(sampleAnnot_2()[,input$SelectVariable2])){
          if (input$ShowDrivers){
            coinertia_plot <- ggplot(data = all_sample_info) +
              geom_point(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], col = colnames(all_sample_info)[9]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], col = colnames(all_sample_info)[9]), shape = 16, size = 3) +
              geom_point(aes_string(x = colnames(all_sample_info)[5], y = colnames(all_sample_info)[6], col = colnames(all_sample_info)[9]), shape = 17, size = 3) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[3], yend = colnames(all_sample_info)[4])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              geom_label_repel(data = axis1_drivers,
                               aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
                               size = 4, segment.size = 0.5,
                               label.padding = unit(0.1, "lines"), label.size = 0.5) +
              theme_bw() +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[1] / sum(mco[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[2] / sum(mco[["mcoa"]]$pseudoeig), 2)))
          }else{

            coinertia_plot <- ggplot(data = all_sample_info) +
              geom_point(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], col = colnames(all_sample_info)[9]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], col = colnames(all_sample_info)[9]), shape = 16, size = 3) +
              geom_point(aes_string(x = colnames(all_sample_info)[5], y = colnames(all_sample_info)[6], col = colnames(all_sample_info)[9]), shape = 17, size = 3) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[3], yend = colnames(all_sample_info)[4])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              # geom_label_repel(data = axis1_drivers,
              #                  aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
              #                  size = 4, segment.size = 0.5,
              #                  label.padding = unit(0.1, "lines"), label.size = 0.5) +
              theme_bw() +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[1] / sum(mco[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[2] / sum(mco[["mcoa"]]$pseudoeig), 2)))
          }

        }else{
          if (input$ShowDrivers){
            coinertia_plot <- ggplot(data = all_sample_info) +
              geom_point(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], col = colnames(all_sample_info)[9]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], col = colnames(all_sample_info)[9]), shape = 16, size = 3) +
              geom_point(aes_string(x = colnames(all_sample_info)[5], y = colnames(all_sample_info)[6], col = colnames(all_sample_info)[9]), shape = 17, size = 3) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[3], yend = colnames(all_sample_info)[4])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              geom_label_repel(data = axis1_drivers,
                               aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
                               size = 4, segment.size = 0.5,
                               label.padding = unit(0.1, "lines"), label.size = 0.5) +
              theme_bw() +
              scale_colour_gradientn(colours = terrain.colors(10)) +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[1] / sum(mco[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[2] / sum(mco[["mcoa"]]$pseudoeig), 2)))
          }else{

            coinertia_plot <- ggplot(data = all_sample_info) +
              geom_point(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], col = colnames(all_sample_info)[9]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], col = colnames(all_sample_info)[9]), shape = 16, size = 3) +
              geom_point(aes_string(x = colnames(all_sample_info)[5], y = colnames(all_sample_info)[6], col = colnames(all_sample_info)[9]), shape = 17, size = 3) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[3], yend = colnames(all_sample_info)[4])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              # geom_label_repel(data = axis1_drivers,
              #                  aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
              #                  size = 4, segment.size = 0.5,
              #                  label.padding = unit(0.1, "lines"), label.size = 0.5) +
              theme_bw() +
              scale_colour_gradientn(colours = terrain.colors(10)) +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[1] / sum(mco[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[2] / sum(mco[["mcoa"]]$pseudoeig), 2)))
          }

        }

        print(plot_grid(coinertia_plot,LegendDF(), ncol = 2, rel_widths = c(1, .3)))

        dev.off()

        fs <- c(fs, paste("Co_inertia_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5, "_df3.",input$selectModule3_p5 , ".pdf", sep = ""))

        # if (!((input$selectModule12_p5 == "Individual_Variables" && is.null(input$selectVariablesExprDat) || (input$selectModule21_p5 == "Individual_Variables" && is.null(input$selectVariables)) ))){
        #   pdf(paste("Heatmap_df1.", input$selectModule12_p5 ,"_df2.", input$selectModule21_p5,".pdf", sep = ""), 10, 10)
        #   heatmap(corr_expr(), margins = c(15,15))
        #   dev.off()
        #   fs <- c(fs, paste("Heatmap_df1.", input$selectModule12_p5 ,"_df2.", input$selectModule21_p5,".pdf", sep = ""))
        # }
        #
        # if (!((input$selectModule13_p5 == "Individual_Variables" && is.null(input$selectVariablesExprDat13) || (input$selectModule31_p5 == "Individual_Variables" && is.null(input$selectVariables31)) ))){
        #   pdf(paste("Heatmap_df1.", input$selectModule13_p5 ,"_df3.", input$selectModule31_p5,".pdf", sep = ""), 10, 10)
        #   heatmap(corr_expr13(), margins = c(15,15))
        #   dev.off()
        #   fs <- c(fs, paste("Heatmap_df1.", input$selectModule13_p5 ,"_df3.", input$selectModule31_p5,".pdf", sep = ""))
        # }
        #
        # if (!((input$selectModule23_p5 == "Individual_Variables" && is.null(input$selectVariablesExprDat23) || (input$selectModule32_p5 == "Individual_Variables" && is.null(input$selectVariables32)) ))){
        #   pdf(paste("Heatmap_df2.", input$selectModule23_p5 ,"_df3.", input$selectModule32_p5,".pdf", sep = ""), 10, 10)
        #   heatmap(corr_expr23(), margins = c(15,15))
        #   dev.off()
        #   fs <- c(fs, paste("Heatmap_df2.", input$selectModule23_p5 ,"_df3.", input$selectModule32_p5,".pdf", sep = ""))
        # }

      }else{
        #
        # svg(paste("HIVEPLOT_", input$traitHive,".svg", sep = ""), width = 10, height = 10)
        #
        # myHive <- hive_data()
        # myHive$nodes$lab <- as.character(myHive$nodes$lab)
        # myHive$nodes$color <- as.character(myHive$nodes$color)
        # myHive$edges$color <- as.character(myHive$edges$color)
        # plotMyHive(myHive, bkgnd = "white", axLabs = c(paste(input$CountingT1, "D1", sep = "_"), paste( input$OmicTable, "D2", sep = "_"), paste(input$OmicTable3, "D3", sep = "_")), axLab.gpar = gpar(col = "#636363"))
        # dev.off()
        #
        # fs <- c(fs,paste("HIVEPLOT_", input$traitHive,".svg", sep = ""))

        svg(paste("Co_inertia_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5, "_df3.",input$selectModule3_p5 , ".svg", sep = ""), width = 10, height = 10)
        all_sample_info <- selected_sampleInfo_all()
        mco <- mcoia()
        axis1_drivers <- selected_axis1_drivers()
        if (!is.numeric(sampleAnnot_2()[,input$SelectVariable2])){
          if (input$ShowDrivers){
            coinertia_plot <- ggplot(data = all_sample_info) +
              geom_point(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], col = colnames(all_sample_info)[9]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], col = colnames(all_sample_info)[9]), shape = 16, size = 3) +
              geom_point(aes_string(x = colnames(all_sample_info)[5], y = colnames(all_sample_info)[6], col = colnames(all_sample_info)[9]), shape = 17, size = 3) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[3], yend = colnames(all_sample_info)[4])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              geom_label_repel(data = axis1_drivers,
                               aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
                               size = 4, segment.size = 0.5,
                               label.padding = unit(0.1, "lines"), label.size = 0.5) +
              theme_bw() +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[1] / sum(mco[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[2] / sum(mco[["mcoa"]]$pseudoeig), 2)))
          }else{
            coinertia_plot <- ggplot(data = all_sample_info) +
              geom_point(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], col = colnames(all_sample_info)[9]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], col = colnames(all_sample_info)[9]), shape = 16, size = 3) +
              geom_point(aes_string(x = colnames(all_sample_info)[5], y = colnames(all_sample_info)[6], col = colnames(all_sample_info)[9]), shape = 17, size = 3) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[3], yend = colnames(all_sample_info)[4])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              # geom_label_repel(data = axis1_drivers,
              #                  aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
              #                  size = 4, segment.size = 0.5,
              #                  label.padding = unit(0.1, "lines"), label.size = 0.5) +
              theme_bw() +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[1] / sum(mco[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[2] / sum(mco[["mcoa"]]$pseudoeig), 2)))
          }

        }else{
          if (input$ShowDrivers){
            coinertia_plot <- ggplot(data = all_sample_info) +
              geom_point(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], col = colnames(all_sample_info)[9]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], col = colnames(all_sample_info)[9]), shape = 16, size = 3) +
              geom_point(aes_string(x = colnames(all_sample_info)[5], y = colnames(all_sample_info)[6], col = colnames(all_sample_info)[9]), shape = 17, size = 3) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[3], yend = colnames(all_sample_info)[4])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              geom_label_repel(data = axis1_drivers,
                               aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
                               size = 4, segment.size = 0.5,
                               label.padding = unit(0.1, "lines"), label.size = 0.5) +
              theme_bw() +
              scale_colour_gradientn(colours = terrain.colors(10)) +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[1] / sum(mco[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[2] / sum(mco[["mcoa"]]$pseudoeig), 2)))
          }else{

            coinertia_plot <- ggplot(data = all_sample_info) +
              geom_point(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], col = colnames(all_sample_info)[9]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], col = colnames(all_sample_info)[9]), shape = 16, size = 3) +
              geom_point(aes_string(x = colnames(all_sample_info)[5], y = colnames(all_sample_info)[6], col = colnames(all_sample_info)[9]), shape = 17, size = 3) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[3], yend = colnames(all_sample_info)[4])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[1], y = colnames(all_sample_info)[2], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              geom_segment(aes_string(x = colnames(all_sample_info)[3], y = colnames(all_sample_info)[4], colour = colnames(all_sample_info)[9], xend = colnames(all_sample_info)[5], yend = colnames(all_sample_info)[6])) +
              # geom_label_repel(data = axis1_drivers,
              #                  aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
              #                  size = 4, segment.size = 0.5,
              #                  label.padding = unit(0.1, "lines"), label.size = 0.5) +
              theme_bw() +
              scale_colour_gradientn(colours = terrain.colors(10)) +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[1] / sum(mco[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mco[["mcoa"]]$pseudoeig[2] / sum(mco[["mcoa"]]$pseudoeig), 2)))
          }

        }

        print(plot_grid(coinertia_plot,LegendDF(), ncol = 2, rel_widths = c(1, .3)))

        dev.off()
        fs <- c(fs, paste("Co_inertia_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5, "_df3.",input$selectModule3_p5 , ".svg", sep = ""))

        # if (!((input$selectModule12_p5 == "Individual_Variables" && is.null(input$selectVariablesExprDat) || (input$selectModule21_p5 == "Individual_Variables" && is.null(input$selectVariables)) ))){
        #   svg(paste("Heatmap_df1.", input$selectModule12_p5 ,"_df2.", input$selectModule21_p5,".svg", sep = ""), 10, 10)
        #   heatmap(corr_expr(), margins = c(15,15))
        #   dev.off()
        #   fs <- c(fs, paste("Heatmap_df1.", input$selectModule12_p5 ,"_df2.", input$selectModule21_p5,".svg", sep = ""))
        # }
        #
        #
        # if (!((input$selectModule13_p5 == "Individual_Variables" && is.null(input$selectVariablesExprDat13) || (input$selectModule31_p5 == "Individual_Variables" && is.null(input$selectVariables31)) ))){
        #   svg(paste("Heatmap_df1.", input$selectModule13_p5 ,"_df3.", input$selectModule31_p5,".svg", sep = ""), 10, 10)
        #   heatmap(corr_expr13(), margins = c(15,15))
        #   dev.off()
        #   fs <- c(fs, paste("Heatmap_df1.", input$selectModule13_p5 ,"_df3.", input$selectModule31_p5,".svg", sep = ""))
        # }
        #
        # if (!((input$selectModule23_p5 == "Individual_Variables" && is.null(input$selectVariablesExprDat23) || (input$selectModule32_p5 == "Individual_Variables" && is.null(input$selectVariables32)) ))){
        #   svg(paste("Heatmap_df2.", input$selectModule23_p5 ,"_df3.", input$selectModule32_p5,".svg", sep = ""), 10, 10)
        #   heatmap(corr_expr23(), margins = c(15,15))
        #   dev.off()
        #   fs <- c(fs, paste("Heatmap_df2.", input$selectModule23_p5 ,"_df3.", input$selectModule32_p5,".svg", sep = ""))
        # }
      }
    }else{
      #-----------------------------------------------------------------------------------------------------

      if (input$pdf_or_svg_p5 == "pdf"){
        # pdf(paste("HIVEPLOT_", input$traitHive,".pdf", sep = ""), width = 10, height = 10)
        #
        # myHive <- hive_data()
        # myHive$nodes$lab <- as.character(myHive$nodes$lab)
        # myHive$nodes$color <- as.character(myHive$nodes$color)
        # myHive$edges$color <- as.character(myHive$edges$color)
        # plotMyHive(myHive, ch = 2, bkgnd = "white", axLabs = c(paste(input$CountingT1, "D1", sep = "_"), paste( input$OmicTable, "D2", sep = "_")), axLab.gpar = gpar(col = "#636363"), rot = c(-90, 90))
        # dev.off()
        #
        # fs <- c(fs,paste("HIVEPLOT_", input$traitHive,".pdf", sep = ""))

        pdf(paste("Co_inertia_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".pdf", sep = ""), width = 10, height = 10)
        if (!is.numeric(sampleAnnot_2()[,input$SelectVariable2])){

          if (input$ShowDrivers){
            coinertia_plot <-
              ggplot(data = selected_sampleInfo_all()) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[7]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[3], y = colnames(selected_sampleInfo_all())[4], col = colnames(selected_sampleInfo_all())[7]), shape = 16, size = 3) +
              scale_shape(solid = TRUE) +
              theme_bw() +
              geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[7], xend = colnames(selected_sampleInfo_all())[3], yend = colnames(selected_sampleInfo_all())[4])) +
              geom_label_repel(data = selected_axis1_drivers(),
                               aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
                               size = 4, segment.size = 0.5,
                               label.padding = unit(0.1, "lines"), label.size = 0.5) +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)))
          }else{
            coinertia_plot <-
              ggplot(data = selected_sampleInfo_all()) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[7]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[3], y = colnames(selected_sampleInfo_all())[4], col = colnames(selected_sampleInfo_all())[7]), shape = 16, size = 3) +
              scale_shape(solid = TRUE) +
              theme_bw() +
              geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[7], xend = colnames(selected_sampleInfo_all())[3], yend = colnames(selected_sampleInfo_all())[4])) +
              # geom_label_repel(data = selected_axis1_drivers(),
              #                  aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
              #                  size = 4, segment.size = 0.5,
              #                  label.padding = unit(0.1, "lines"), label.size = 0.5) +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)))
          }
        }else{
          if (input$ShowDrivers){
            coinertia_plot <-
              ggplot(data = selected_sampleInfo_all()) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[7]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[3], y = colnames(selected_sampleInfo_all())[4], col = colnames(selected_sampleInfo_all())[7]), shape = 16, size = 3) +
              scale_shape(solid = TRUE) +
              theme_bw() +
              scale_colour_gradientn(colours = terrain.colors(10)) +
              geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[7], xend = colnames(selected_sampleInfo_all())[3], yend = colnames(selected_sampleInfo_all())[4])) +
              geom_label_repel(data = selected_axis1_drivers(),
                               aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
                               size = 4, segment.size = 0.5,
                               label.padding = unit(0.1, "lines"), label.size = 0.5) +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)))
          }else{
            coinertia_plot <-
              ggplot(data = selected_sampleInfo_all()) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[7]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[3], y = colnames(selected_sampleInfo_all())[4], col = colnames(selected_sampleInfo_all())[7]), shape = 16, size = 3) +
              scale_shape(solid = TRUE) +
              theme_bw() +
              scale_colour_gradientn(colours = terrain.colors(10)) +
              geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[7], xend = colnames(selected_sampleInfo_all())[3], yend = colnames(selected_sampleInfo_all())[4])) +
              # geom_label_repel(data = selected_axis1_drivers(),
              #                  aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
              #                  size = 4, segment.size = 0.5,
              #                  label.padding = unit(0.1, "lines"), label.size = 0.5) +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)))
          }

        }

        bivariate_plot <-
          ggplot(data = bivariate_plot(),
                 aes_string(x = colnames(bivariate_plot())[1],
                            y = colnames(bivariate_plot())[2],
                            label = "sampleName")) +
          theme_bw() +
          geom_point() +
          geom_abline(intercept = 0, slope = 1, col = "red") +
          geom_vline(xintercept = 0) +
          geom_hline(yintercept = 0) +
          geom_text(vjust = 0, nudge_y = 0.02) +
          annotate("text", label = "maximal correlation", x = -0.5, y =-0.5, color = "red", size = 5)

        print(plot_grid(coinertia_plot,LegendDF(),bivariate_plot, ncol = 2, rel_widths = c(1, .3)))
        #
        # print(ggdendrogram(distTree_coia()))
        #
        # print(ggplot(data = distance_mcoia_df_supp_info(), aes(x= sampleName, y = distance_coia.v, col = color)) +
        #         geom_point() +
        #         theme(axis.text.x = element_text(angle = 90, hjust = 1)))
        #
        # print(ggplot(distance_mcoia_df_supp_info(), aes(distance_coia.v, fill = color)) + geom_histogram(binwidth = 0.1))

        dev.off()

        fs <- c(fs, paste("Co_inertia_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".pdf", sep = ""))



        pdf(paste("Procrustes_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".pdf", sep = ""), width = 10, height = 10)
        print(ggplot(data = sample_info_PA()) + geom_point(aes_string(x = 'Axis1', y= 'Axis2', col = input$SelectVariable2), size = 3) +
                geom_point(aes_string(x = 'Axis.1', y= 'Axis.2', col = input$SelectVariable2), size = 3, shape = 5) +
                geom_segment(aes_string(x = 'Axis1', y = 'Axis2', xend = 'Axis.1', yend = 'Axis.2', colour = input$SelectVariable2)))

        dev.off()
        fs <- c(fs, paste("Procrustes_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".pdf", sep = ""))

        # if ((input$selectModule12_p5 == "Individual_Variables" && is.null(input$selectVariablesExprDat) || (input$selectModule21_p5 == "Individual_Variables" && is.null(input$selectVariables)) )){
        #   print("something")
        #
        # }else{
        #   pdf(paste("Heatmap_df1.", input$selectModule12_p5 ,"_df2.", input$selectModule21_p5,".pdf", sep = ""), 10, 10)
        #   heatmap(corr_expr(), ColSideColors = my_col, margins = c(15,15))
        #   dev.off()
        #   fs <- c(fs, paste("Heatmap_df1.", input$selectModule12_p5 ,"_df2.", input$selectModule21_p5,".pdf", sep = ""))
        # }

      }else{
        # svg(paste("HIVEPLOT_", input$traitHive,".svg", sep = ""), width = 10, height = 10)
        #
        # myHive <- hive_data()
        # myHive$nodes$lab <- as.character(myHive$nodes$lab)
        # myHive$nodes$color <- as.character(myHive$nodes$color)
        # myHive$edges$color <- as.character(myHive$edges$color)
        # plotMyHive(myHive, ch = 2, bkgnd = "white", axLabs = c(paste(input$CountingT1, "D1", sep = "_"), paste( input$OmicTable, "D2", sep = "_")), axLab.gpar = gpar(col = "#636363"), rot = c(-90, 90))
        # dev.off()
        #
        # fs <- c(fs,paste("HIVEPLOT_", input$traitHive,".svg", sep = ""))
        svg(paste("Co_inertia_mainPlot_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".svg", sep = ""), width = 10, height = 10)
        if (!is.numeric(sampleAnnot_2()[,input$SelectVariable2])){
          if (input$ShowDrivers){
            coinertia_plot <-
              ggplot(data = selected_sampleInfo_all()) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[7]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[3], y = colnames(selected_sampleInfo_all())[4], col = colnames(selected_sampleInfo_all())[7]), shape = 16, size = 3) +
              scale_shape(solid = TRUE) +
              theme_bw() +
              geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[7], xend = colnames(selected_sampleInfo_all())[3], yend = colnames(selected_sampleInfo_all())[4])) +
              geom_label_repel(data = selected_axis1_drivers(),
                               aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
                               size = 4, segment.size = 0.5,
                               label.padding = unit(0.1, "lines"), label.size = 0.5) +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)))
          }else{
            coinertia_plot <-
              ggplot(data = selected_sampleInfo_all()) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[7]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[3], y = colnames(selected_sampleInfo_all())[4], col = colnames(selected_sampleInfo_all())[7]), shape = 16, size = 3) +
              scale_shape(solid = TRUE) +
              theme_bw() +
              geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[7], xend = colnames(selected_sampleInfo_all())[3], yend = colnames(selected_sampleInfo_all())[4])) +
              # geom_label_repel(data = selected_axis1_drivers(),
              #                  aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
              #                  size = 4, segment.size = 0.5,
              #                  label.padding = unit(0.1, "lines"), label.size = 0.5) +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)))
          }
        }else{
          if (input$ShowDrivers){
            coinertia_plot <-
              ggplot(data = selected_sampleInfo_all()) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[7]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[3], y = colnames(selected_sampleInfo_all())[4], col = colnames(selected_sampleInfo_all())[7]), shape = 16, size = 3) +
              scale_shape(solid = TRUE) +
              theme_bw() +
              scale_colour_gradientn(colours = terrain.colors(10)) +
              geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[7], xend = colnames(selected_sampleInfo_all())[3], yend = colnames(selected_sampleInfo_all())[4])) +
              geom_label_repel(data = selected_axis1_drivers(),
                               aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
                               size = 4, segment.size = 0.5,
                               label.padding = unit(0.1, "lines"), label.size = 0.5) +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)))
          }else{
            coinertia_plot <-
              ggplot(data = selected_sampleInfo_all()) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[7]), size = 3, shape = 15) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[3], y = colnames(selected_sampleInfo_all())[4], col = colnames(selected_sampleInfo_all())[7]), shape = 16, size = 3) +
              scale_shape(solid = TRUE) +
              theme_bw() +
              scale_colour_gradientn(colours = terrain.colors(10)) +
              geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[7], xend = colnames(selected_sampleInfo_all())[3], yend = colnames(selected_sampleInfo_all())[4])) +
              # geom_label_repel(data = selected_axis1_drivers(),
              #                  aes(x = 0.10 * SV1, y = 0.10 * SV2, label = spec, fill = feature),
              #                  size = 4, segment.size = 0.5,
              #                  label.padding = unit(0.1, "lines"), label.size = 0.5) +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)))
          }

        }

        bivariate_plot <-
          ggplot(data = bivariate_plot(),
                 aes_string(x = colnames(bivariate_plot())[1],
                            y = colnames(bivariate_plot())[2],
                            label = "sampleName")) +
          theme_bw() +
          geom_point() +
          geom_abline(intercept = 0, slope = 1, col = "red") +
          geom_vline(xintercept = 0) +
          geom_hline(yintercept = 0) +
          geom_text(vjust = 0, nudge_y = 0.02) +
          annotate("text", label = "maximal correlation", x = -0.5, y =-0.5, color = "red", size = 5)
        print(plot_grid(coinertia_plot,LegendDF(),bivariate_plot, ncol = 2, rel_widths = c(1, .3)))
        dev.off()
        fs <- c(fs, paste("Co_inertia_mainPlot_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".svg", sep = ""))
        svg(paste("Co_inertia_bivariatePlot_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".svg", sep = ""), width = 10, height = 10)
        print(ggplot(data = bivariate_plot(),
                     aes_string(x = colnames(bivariate_plot())[1],
                                y = colnames(bivariate_plot())[2],
                                label = "sampleName")) +
                geom_point() +
                geom_abline(intercept = 0, slope = 1, col = "red") +
                geom_vline(xintercept = 0) +
                geom_hline(yintercept = 0) +
                geom_text(vjust = 0, nudge_y = 0.02) +
                annotate("text", label = "maximal correlation", x = -0.5, y =-0.5, color = "red", size = 5))
        dev.off()
        fs <- c(fs, paste("Co_inertia_bivariatePlot_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".svg", sep = ""))

        svg(paste("Co_inertia_dendrogramme_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".svg", sep = ""), width = 10, height = 10)
        print(ggdendrogram(distTree_coia()))
        dev.off()
        fs <- c(fs, paste("Co_inertia_dendrogramme_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".svg", sep = ""))

        svg(paste("Co_inertia_distPoint_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".svg", sep = ""), width = 10, height = 10)
        print(ggplot(data = distance_mcoia_df_supp_info(), aes(x= sampleName, y = distance_coia.v, col = color)) +
                geom_point() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)))
        dev.off()
        fs <- c(fs, paste("Co_inertia_distPoint_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".svg", sep = ""))


        svg(paste("Co_inertia_distHist_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".svg", sep = ""), width = 10, height = 10)
        print(ggplot(distance_mcoia_df_supp_info(), aes(distance_coia.v, fill = color)) + geom_histogram(binwidth = 0.1))

        dev.off()

        fs <- c(fs, paste("Co_inertia_distHist_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".svg", sep = ""))



        svg(paste("Procrustes_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".svg", sep = ""), width = 10, height = 10)
        print(ggplot(data = sample_info_PA()) + geom_point(aes_string(x = 'Axis1', y= 'Axis2', col = input$SelectVariable2), size = 3) +
                geom_point(aes_string(x = 'Axis.1', y= 'Axis.2', col = input$SelectVariable2), size = 3, shape = 5) +
                geom_segment(aes_string(x = 'Axis1', y = 'Axis2', xend = 'Axis.1', yend = 'Axis.2', colour = input$SelectVariable2)))

        dev.off()
        fs <- c(fs, paste("Procrustes_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".svg", sep = ""))

        # if ((input$selectModule12_p5 == "Individual_Variables" && is.null(input$selectVariablesExprDat) || (input$selectModule21_p5 == "Individual_Variables" && is.null(input$selectVariables)) )){
        #
        #
        # }else{
        #   svg(paste("Heatmap_df1.", input$selectModule12_p5 ,"_df2.", input$selectModule21_p5,".svg", sep = ""), 10, 10)
        #   heatmap(corr_expr(), ColSideColors = my_col, margins = c(15,15))
        #   dev.off()
        #   fs <- c(fs, paste("Heatmap_df1.", input$selectModule12_p5 ,"_df2.", input$selectModule21_p5,".svg", sep = ""))
        # }

      }
    }


    zip(zipfile=filename, files=fs)


  },
  contentType = "application/zip")

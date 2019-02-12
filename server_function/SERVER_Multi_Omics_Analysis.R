#################################################
#### SERVER FUNCTIONS : MULTI-OMICS ANALYSIS ####
#################################################
# This script gathers the functions related to the multi-omics analysis tabulation of MiBiOmics. It contains all the
# functions needed for the co-inertia analysis, procrustes analysis and the association of modules extracted from
# WGCNA.

#### REACTIVE OBJECTS ####

mcoia <- reactive({

  selectMetab <- selectedMetab()
  # no_variation <- 0
  #
  # for (i in 1:nrow(selected_moduleOTU())){
  #   if (length(unique(selected_moduleOTU()[i,])) == 1){
  #     no_variation <- no_variation + 1
  #   }
  # }
  # for (i in 1:nrow(selectedMetab())){
  #   if (length(unique(selectedMetab()[i,])) == 1){
  #     no_variation <- no_variation + 1
  #   }
  # }
  #
  # if (no_variation == 0){
  df1_df2 <- list()
  df1_df2[["df1"]] <- t(as.data.frame(selected_moduleOTU()))
  df1_df2[["df2"]] <- t(as.data.frame(selectMetab))
  mcia(df1_df2, nsc = F)
  # }

})


selected_sampleInfo <- reactive({
  conditionCol <- input$SelectVariable2
  vector_info1 <- c(rep(sampleAnnot_3()[,conditionCol], 2))
  sample_info_mcoia <- data.frame(as.data.frame(mcoia()[["mcoa"]]$Tl1), as.data.frame(vector_info1))
  sample_info_mcoia
})

selected_sampleInfo_df1 <- reactive({
  sample_info_mcoia_df1 <-selected_sampleInfo()[1:(nrow(selected_sampleInfo())/2),]
  sample_info_mcoia_df1
})


selected_sampleInfo_df2 <- reactive({
  sample_info_mcoia_df2 <-selected_sampleInfo()[(nrow(selected_sampleInfo())/2+1):(nrow(selected_sampleInfo())),]
  sample_info_mcoia_df2
})


selected_sampleInfo_all <- reactive({
  sample_info_mcoia_all <- data.frame(selected_sampleInfo_df1(), selected_sampleInfo_df2())
  colnames(sample_info_mcoia_all) <- c("col1", "col2", "col3", "col4", "col5", "col6")
  sample_info_mcoia_all
})

### A dataframe containing the information about the length of each vector in the co-inertia
distance_mcoia_df <- reactive({
  distance_coia.v <- c()

  for (i in 1:nrow(selected_sampleInfo_all())){
    # Calculate distance between two point in co-inertia
    dist_coia <- sqrt((selected_sampleInfo_all()[i, 4] - selected_sampleInfo_all()[i, 1])^2 + (selected_sampleInfo_all()[i, 5] - selected_sampleInfo_all()[i, 2])^2)
    distance_coia.v <- c(distance_coia.v, dist_coia)
  }

  distance_mcoia.df <- data.frame(distance_coia.v)
  rownames(distance_mcoia.df) <- rownames(sampleAnnot_3())
  distance_mcoia.df$sampleName <- rownames(sampleAnnot_3())
  distance_mcoia.df

})

distTree_coia <- reactive({
  hclust(dist(distance_mcoia_df()), method = "average")
})

distance_mcoia_df_supp_info <- reactive({
  distance_mcoia.df <- distance_mcoia_df()
  distance_mcoia.df$color <- sampleAnnot_3()[,input$SelectVariable2]

  distance_mcoia.df <- distance_mcoia.df[order(distance_mcoia.df$distance_coia.v),]
  distance_mcoia.df$sampleName <- as.character(distance_mcoia.df$sampleName)
  distance_mcoia.df$sampleName <- factor(distance_mcoia.df$sampleName, levels= unique(distance_mcoia.df$sampleName))
  distance_mcoia.df
})

bivariate_plot <- reactive({
  bivariate_plot <- data.frame(mcoia()[["mcoa"]]$Tl1[1:(nrow(selected_sampleInfo())/2), 1], mcoia()[["mcoa"]]$Tl1[(nrow(selected_sampleInfo())/2+1):(nrow(selected_sampleInfo())), 1])
  if (input$LoadExample2 == "Yes"){
    colnames(bivariate_plot) <- c("Coinertia_OTUs_axis", "Coinertia_Metabolite_axis")
  }else{
    colnames(bivariate_plot) <- c(paste("Coinertia_", input$CountingT1, "_axis", sep =""), paste("Coinertia_", input$OmicTable, "D2_axis", sep=""))
  }
  bivariate_plot$sampleName <- rownames(sampleAnnot_3())
  bivariate_plot
})

RV_score <- reactive({
  mcoia()[["mcoa"]]$RV[1,2]
})

selected_axis1_drivers <- reactive({
  axis1_drivers <- data.frame(Variables_PosX = character(), Variables_PosY = character())
  for (i in 1:nrow(mcoia()[["mcoa"]]$Tco)){
    if (mcoia()[["mcoa"]]$Tco$SV1[i] < quantile(mcoia()[["mcoa"]]$Tco$SV1, 0.01)){
      axis1_drivers <- rbind(axis1_drivers, mcoia()[["mcoa"]]$Tco[i,])
    }
    if (mcoia()[["mcoa"]]$Tco$SV1[i] > quantile(mcoia()[["mcoa"]]$Tco$SV1, 0.99)){
      axis1_drivers <- rbind(axis1_drivers, mcoia()[["mcoa"]]$Tco[i,])
    }
  }
  print(axis1_drivers)
  OTUs <- c()
  for (i in 1:nrow(axis1_drivers)){
    OTUs <- c(OTUs, substr(rownames(axis1_drivers)[i], 1, nchar(rownames(axis1_drivers)[i])-4))
  }
  axis1_drivers$OTUs <- OTUs
  feature_type_axis1 <- c()
  # If taxa file <- change OTUs sequence to family name
  if (input$TaxonFile1 || input$LoadExample2 == "Yes"){
    if (input$OmicTable == "OTUs" && input$CountingT1 == "OTUs"){
      taxon_interest_axis1 <- data.frame(taxTable2()[which(taxTable2()[["rn"]] %in% axis1_drivers$OTUs),])
    }else{
      taxon_interest_axis1 <- data.frame(taxTable2()[which(rownames(taxTable2()) %in% rownames(axis1_drivers)),])
    }

    j <- 1
  }
  for (i in 1:nrow(axis1_drivers)){
    if (rownames(axis1_drivers)[i] %in% c(colnames(exprDat_3()),colnames(exprDatSec_2()))){
      if(input$LoadExample2 == "Yes"){
        feature_type_axis1 <- c(feature_type_axis1, "OTUs")
        axis1_drivers$Spec[i] <- taxon_interest_axis1[[input$selectTaxonomy]][j]
        j <- j + 1
      }else{

        if (input$TaxonFile1 && (rownames(axis1_drivers)[i] %in% colnames(exprDat_3()))){
          print("inhere")
          feature_type_axis1 <- c(feature_type_axis1, input$CountingT1)
          axis1_drivers$Spec[i] <- taxon_interest_axis1[[input$selectTaxonomy]][j]
          j <- j + 1
        }else{
          feature_type_axis1 <- c(feature_type_axis1, input$OmicTable)
          axis1_drivers$Spec[i] <- rownames(axis1_drivers)[i]
        }
      }
    }else{
      if(input$LoadExample2 == "Yes"){
        axis1_drivers$Spec[i] <- rownames(axis1_drivers)[i]
        feature_type_axis1 <- c(feature_type_axis1, "metabolites")
      }else{
        if (input$CountingT1 == "OTUs" && input$OmicTable == "OTUs" && input$LoadExample2 == "No"){
          feature_type_axis1 <- rep("OTUs", nrow(axis1_drivers))
          Spec <- c()
          for (i in 1:nrow(axis1_drivers)){
            Spec <- c(Spec, taxon_interest_axis1[which(taxon_interest_axis1$rn == axis1_drivers$OTUs[i]), input$selectTaxonomy])

          }
          axis1_drivers$Spec <- Spec
        }
        #axis1_drivers$Spec[i] <- rownames(axis1_drivers)[i]
        #feature_type_axis1 <- c(feature_type_axis1, input$OmicTable)
      }
    }

  }
  axis1_drivers$feature_type_axis1 <- feature_type_axis1
  axis1_drivers
})


#### PROCRUSTES


dudi_exprDat <- reactive({

  if (input$selectModule1_p5 == "General"){
    dudi_exprDat <- dudi.pca(exprDat_3(), scannf = FALSE, nf = 2 )
  }else{
    Module_OTU <- exprDat_3()[,modGenes]
    dudi_exprDat <- dudi.pca(Module_OTU, scannf = FALSE, nf = 2 )
  }

  dudi_exprDat
})


dudi_exprDatSec <- reactive({

  dudi_exprDatSec <- dudi.pca(selectedMetab(), scannf = FALSE, nf = 2 )

  dudi_exprDatSec
})

PA_object <- reactive({
  PA_object <- procrustes(dudi_exprDat(), dudi_exprDatSec(), scale = TRUE)
  PA_object
})


sample_info_PA <- reactive({
  sample_info_X_PA <- data.frame(PA_object()[["X"]], sampleAnnot_3())
  colnames(sample_info_X_PA) <- c("Axis1", "Axis2", colnames(sampleAnnot_3()))
  sample_info_Y_PA <- data.frame(PA_object()[["Yrot"]], sampleAnnot_3())
  colnames(sample_info_Y_PA) <- c("Axis.1", "Axis.2", colnames(sampleAnnot_3()))
  sample_info_PA <- data.frame(sample_info_X_PA, sample_info_Y_PA)
  sample_info_PA
})

#### PRE - HEATMAP

corr_MEs_D1D2 <- reactive({
  corr_MEs_D1D2 <- cor(selectedMEs(), selectedMEs2(), use = "p", method = "spearman")
  corr_MEs_D1D2 <- as.data.frame(corr_MEs_D1D2)
  row.order <- hclust(dist(corr_MEs_D1D2, method = "euclidean"), method = "ward.D")$order
  col.order <- hclust(dist(t(corr_MEs_D1D2), method = "euclidean"), method = "ward.D")$order
  corr_MEs_D1D2 <- corr_MEs_D1D2[row.order, col.order]
  t(corr_MEs_D1D2)
})

p_val_MEs_D1D2 <- reactive({
  pvalue_MEs_D1_D2 <- corPvalueStudent(cor(selectedMEs(), selectedMEs2(), use = "p", method = "spearman"), nrow(selectedMEs()))
  pvalue_MEs_D1_D2 <- as.data.frame(pvalue_MEs_D1_D2)
  row.order <- hclust(dist(pvalue_MEs_D1_D2, method = "euclidean"), method = "ward.D")$order
  col.order <- hclust(dist(t(pvalue_MEs_D1_D2), method = "euclidean"), method = "ward.D")$order
  pvalue_MEs_D1_D2 <- pvalue_MEs_D1_D2[row.order, col.order]
  t(pvalue_MEs_D1_D2)
})

corr_expr <- reactive({
  # no_variation <- 0
  #
  # for (i in 1:nrow(exprDat_4())){
  #   if (length(unique(exprDat_4()[i,])) == 1){
  #     no_variation <- no_variation + 1
  #   }
  # }
  # for (i in 1:nrow(exprDatSec_4())){
  #   if (length(unique(exprDatSec_4()[i,])) == 1){
  #     no_variation <- no_variation + 1
  #   }
  # }
  # validate(
  #   need(no_variation == 0, "Blablabla")
  # )
  corr_expr_1 <- cor(exprDat_4(), exprDatSec_5(), use = "p", method = "spearman")
  corr_expr_1 <- as.data.frame(corr_expr_1)
  row.order <- hclust(dist(corr_expr_1, method = "euclidean"), method = "ward.D")$order
  col.order <- hclust(dist(t(corr_expr_1), method = "euclidean"), method = "ward.D")$order
  corr_expr_new <- corr_expr_1[row.order, col.order]
  # taxon_annotation <- taxTable1()[which(taxTable1()[["rn"]] %in% rownames(corr_expr_new))]
  # taxon_group <- as.vector(taxon_annotation[["Family"]])
  # if (nchar(rownames(corr_expr_new)[5]) >20){
  #   for (row in 1:nrow(corr_expr_new)){
  #     rownames(corr_expr_new)[row] <- paste("OTU", row, "_", taxon_group[row], sep="")
  #   }
  # }
  t(corr_expr_new)
})

family_group <- reactive({
  taxon_group <- c()
  for (i in 1:ncol(corr_expr())){
    taxon_group <- c(taxon_group, (taxTable2()[which(taxTable2()[["rn"]] == colnames(corr_expr())[i]),input$selectTaxonomy]))
  }
  taxon_group
})

family_group_D2 <- reactive({
  if (input$OmicTable == "OTUs"){
    taxon_group <- c()
    for (i in 1:nrow(corr_expr())){
      print(rownames(corr_expr())[i])
      print(taxTable2()[which(taxTable2()[["rn"]] == rownames(corr_expr())[i]),input$selectTaxonomy])
      taxon_group <- c(taxon_group, (taxTable2()[which(taxTable2()[["rn"]] == rownames(corr_expr())[i]),input$selectTaxonomy]))
    }
    print(length(taxon_group))
    print(nrow(corr_expr()))
    taxon_group
  }

})

#### BI-PARTITE NETWORK
adjacency_expr <- reactive({
  corr_expr_1 <- cor(exprDat_4(), exprDatSec_5(), use = "p", method = "spearman")
  corr_expr_1 <- as.data.frame(corr_expr_1)
  #corr_expr_1[which(corr_expr_1 < 0.30 && corr_expr_1 > -0.30),] <- 0
  corr_expr_1[abs(corr_expr_1) < 0.5] <- 0
  #corr_expr_1[(corr_expr_1 < 0.70 && corr_expr_1 > 0.50) || (corr_expr_1 > -0.70 && corr_expr_1 < -0.50)] <- 1
  corr_expr_1[abs(corr_expr_1)>=0.5] <- 1
  #corr_expr_1 <- t(corr_expr_1)
  if (nchar(rownames(corr_expr_1)[5]) >20){
    for (row in 1:nrow(corr_expr_1)){

      rownames(corr_expr_1)[row] <- paste(input$CountingT1, "_D1", row, sep="")
    }
  }
  if (nchar(colnames(corr_expr_1)[5]) >20){
    for (col in 1:ncol(corr_expr_1)){

      colnames(corr_expr_1)[col] <- paste(input$OmicTable, "_D2", row, sep="")
    }
  }
  adjacency_expr <- corr_expr_1
  adjacency_expr
})

bip <- reactive({

  bip = network::network(adjacency_expr(),
                         matrix.type = "bipartite",
                         ignore.eval = FALSE,
                         names.eval = "weights")


  bip
})


#### INTERFACE VARIABLES ####


output$SelectModule1 <- renderUI({
  selectInput("selectModule1_p5",
              label= "Choose a module color: ",
              choices = c("General", substr(names(selectedMEs()), 3, 40)),
              selected = "Individual Variables")
})

output$SelectModule11 <- renderUI({
  selectInput("selectModule11_p5",
              label= "Choose a module color: ",
              choices = c("Individual_Variables", substr(names(selectedMEs()), 3, 40)),
              selected = "Individual_Variables")
})

output$SelectModule2 <- renderUI({
  selectInput("selectModule2_p5",
              label= "Choose a module color: ",
              choices = c("General", substr(names(selectedMEs2()), 3, 40)),
              selected = "General")
})

output$SelectModule22 <- renderUI({
  selectInput("selectModule22_p5",
              label= "Choose a module color: ",
              choices = c("Individual_Variables", substr(names(selectedMEs2()), 3, 40)),
              selected = "Individual_Variables")
})

output$SelectVariable2 <- renderUI({
  selectInput("SelectVariable2",
              label = "Choose a variable: ",
              choices = colnames(sampleAnnot_2()),
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
  if (!is.numeric(sampleAnnot_2()[,input$SelectVariable2])){
    if (input$ShowDrivers){
      coinertia_plot <-
        ggplot(data = selected_sampleInfo_all()) +
        geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[3]), size = 3) +
        geom_point(aes_string(x = colnames(selected_sampleInfo_all())[4], y = colnames(selected_sampleInfo_all())[5], col = colnames(selected_sampleInfo_all())[3]), shape = 5, size = 3) +
        scale_shape(solid = TRUE) +
        theme_bw() +
        geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[3], xend = colnames(selected_sampleInfo_all())[4], yend = colnames(selected_sampleInfo_all())[5])) +
        geom_label_repel(data = selected_axis1_drivers(),
                         aes(x = 0.10 * SV1, y = 0.10 * SV2, label = Spec, fill = feature_type_axis1),
                         size = 4, segment.size = 0.5,
                         label.padding = unit(0.1, "lines"), label.size = 0.5) +
        labs(x = sprintf("Axis1 [%s%% Variance]",
                         100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
             y = sprintf("Axis2 [%s%% Variance]",
                         100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)))
    }else{
      coinertia_plot <-
        ggplot(data = selected_sampleInfo_all()) +
        geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[3]), size = 3) +
        geom_point(aes_string(x = colnames(selected_sampleInfo_all())[4], y = colnames(selected_sampleInfo_all())[5], col = colnames(selected_sampleInfo_all())[3]), shape = 5, size = 3) +
        scale_shape(solid = TRUE) +
        theme_bw() +
        geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[3], xend = colnames(selected_sampleInfo_all())[4], yend = colnames(selected_sampleInfo_all())[5])) +
        # geom_label_repel(data = selected_axis1_drivers(),
        #                  aes(x = 0.10 * SV1, y = 0.10 * SV2, label = Spec, fill = feature_type_axis1),
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
        geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[3]), size = 3) +
        geom_point(aes_string(x = colnames(selected_sampleInfo_all())[4], y = colnames(selected_sampleInfo_all())[5], col = colnames(selected_sampleInfo_all())[3]), shape = 5, size = 3) +
        scale_shape(solid = TRUE) +
        geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[3], xend = colnames(selected_sampleInfo_all())[4], yend = colnames(selected_sampleInfo_all())[5])) +
        geom_label_repel(data = selected_axis1_drivers(),
                         aes(x = 0.10 * SV1, y = 0.10 * SV2, label = Spec, fill = feature_type_axis1),
                         size = 4, segment.size = 0.5,
                         label.padding = unit(0.1, "lines"), label.size = 0.5) +
        theme_bw() +
        scale_colour_gradientn(colours = terrain.colors(10)) +
        labs(x = sprintf("Axis1 [%s%% Variance]",
                         100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
             y = sprintf("Axis2 [%s%% Variance]",
                         100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)))
    }else{
      coinertia_plot <-
        ggplot(data = selected_sampleInfo_all()) +
        geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[3]), size = 3) +
        geom_point(aes_string(x = colnames(selected_sampleInfo_all())[4], y = colnames(selected_sampleInfo_all())[5], col = colnames(selected_sampleInfo_all())[3]), shape = 5, size = 3) +
        scale_shape(solid = TRUE) +
        geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[3], xend = colnames(selected_sampleInfo_all())[4], yend = colnames(selected_sampleInfo_all())[5])) +
        # geom_label_repel(data = selected_axis1_drivers(),
        #                  aes(x = 0.10 * SV1, y = 0.10 * SV2, label = Spec, fill = feature_type_axis1),
        #                  size = 4, segment.size = 0.5,
        #                  label.padding = unit(0.1, "lines"), label.size = 0.5) +
        theme_bw() +
        scale_colour_gradientn(colours = terrain.colors(10)) +
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
  plot_grid(bivariate_plot, coinertia_plot, labels = c("A", "B"), nrow = 2, align = "v")
})

output$coinertia_bivariate <- renderPlot({

  dendrogramme_coia <-
    ggdendrogram(distTree_coia())
  plot_dist_coia <-
    ggplot(data = distance_mcoia_df_supp_info(), aes(x= sampleName, y = distance_coia.v, col = color)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot_hist_coia <-
    ggplot(distance_mcoia_df_supp_info(), aes(distance_coia.v, fill = color)) + geom_histogram(binwidth = 0.1)
  plot_grid(dendrogramme_coia,
            plot_dist_coia,
            plot_hist_coia,
            labels = c("A", "B", "C"),
            nrow = 2,
            align = "h")

})

output$RV <- renderPrint({
  RV_score()
})

output$RV2 <- renderPrint({
  RV_score()
})

output$PA <- renderPlot({
  ggplot(data = sample_info_PA()) + geom_point(aes_string(x = 'Axis1', y= 'Axis2', col = input$SelectVariable2), size = 3) +
    geom_point(aes_string(x = 'Axis.1', y= 'Axis.2', col = input$SelectVariable2), size = 3, shape = 5) +
    geom_segment(aes_string(x = 'Axis1', y = 'Axis2', xend = 'Axis.1', yend = 'Axis.2', colour = input$SelectVariable2))
})


output$HEATMAP_MEs <- renderIheatmap({
  main_heatmap(corr_MEs_D1D2(), name = "Correlation between MEs") %>%
    add_row_labels() %>%
    add_col_labels() %>%
    add_col_clustering() %>%
    add_row_clustering() %>%
    add_col_title("Modules of the first dataset") %>%
    add_row_title("Modules of the second dataset")

})

output$HEATMAP <- renderIheatmap({
  # if (input$OmicTable == "OTUs"){
  #   main_heatmap(corr_expr(), name = "Correlation") %>%
  #     add_row_labels() %>%
  #     add_col_clustering() %>%
  #     add_row_clustering() %>%
  #     add_col_title("Variables of the first datasets in each module") %>%
  #     add_row_title("Variables of the second dataset") %>%
  #     add_col_annotation(data.frame("Groups" = family_group())) %>%
  #     add_row_annotation(data.frame("Row Groups" = family_group_D2()))
  # }else{
  main_heatmap(corr_expr(), name = "Correlation") %>%
    add_row_labels() %>%
    add_col_clustering() %>%
    add_row_clustering() %>%
    add_col_title("Variables of the first datasets in each module") %>%
    add_row_title("Variables of the second dataset") %>%
    add_col_annotation(data.frame("Groups" = family_group()))
  # }

})

output$bipartite_network <- renderPlot({
  col = c("actor" = "red", "event" = "blue")
  ggnet2(bip(), color = "mode", palette = col, label = TRUE, layout.exp = 0.25, shape= "mode") +
    theme(legend.position="none")

})

#### DOWNLOADS ####

output$Download_Coinertia_drivers <- downloadHandler(
  filename = function(){
    paste(input$selectModule1_p5, "_", input$selectModule2_p5, "coinertia_drivers.csv", sep="")
  },
  content = function (filename){
    write.csv(as.data.frame(mcoia()[["mcoa"]]$Tco), filename)
  }
)

output$Download_Multivariate_Analysis <- downloadHandler(
  filename = function(){
    paste("MULTIVARIATE_ANALYSIS.zip")
  },
  content = function (filename){
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    if (input$pdf_or_svg_p5 == "pdf"){
      pdf(paste("Co_inertia_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".pdf", sep = ""), width = 10, height = 10)
      print(ggplot(data = selected_sampleInfo_all()) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[3]), size = 3) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[4], y = colnames(selected_sampleInfo_all())[5], col = colnames(selected_sampleInfo_all())[3]), shape = 5, size = 3) +
              scale_shape(solid = TRUE) +
              geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[3], xend = colnames(selected_sampleInfo_all())[4], yend = colnames(selected_sampleInfo_all())[5])) +
              geom_label_repel(data = selected_axis1_drivers(),
                               aes(x = 0.10 * SV1, y = 0.10 * SV2, label = Spec, fill = feature_type_axis1),
                               size = 4, segment.size = 0.5,
                               label.padding = unit(0.1, "lines"), label.size = 0.5) +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2))))
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

      print(ggdendrogram(distTree_coia()))

      print(ggplot(data = distance_mcoia_df_supp_info(), aes(x= sampleName, y = distance_coia.v, col = color)) +
              geom_point() +
              theme(axis.text.x = element_text(angle = 90, hjust = 1)))

      print(ggplot(distance_mcoia_df_supp_info(), aes(distance_coia.v, fill = color)) + geom_histogram(binwidth = 0.1))

      dev.off()

      fs <- c(fs, paste("Co_inertia_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".pdf", sep = ""))



      pdf(paste("Procrustes_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".pdf", sep = ""), width = 10, height = 10)
      print(ggplot(data = sample_info_PA()) + geom_point(aes_string(x = 'Axis1', y= 'Axis2', col = input$SelectVariable2), size = 3) +
              geom_point(aes_string(x = 'Axis.1', y= 'Axis.2', col = input$SelectVariable2), size = 3, shape = 5) +
              geom_segment(aes_string(x = 'Axis1', y = 'Axis2', xend = 'Axis.1', yend = 'Axis.2', colour = input$SelectVariable2)))

      dev.off()
      fs <- c(fs, paste("Procrustes_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".pdf", sep = ""))

      if ((input$selectModule11_p5 == "Individual_Variables" && is.null(input$selectVariablesExprDat) || (input$selectModule22_p5 == "Individual_Variables" && is.null(input$selectVariables)) )){
        print("something")

      }else{
        pdf(paste("Heatmap_df1.", input$selectModule11_p5 ,"_df2.", input$selectModule22_p5,".pdf", sep = ""), 10, 10)
        my_group=as.numeric(as.factor(family_group()))
        my_col=brewer.pal(9, "Set1")[my_group]
        heatmap(corr_expr(), ColSideColors = my_col, margins = c(15,15))

        legend("topright",
               legend = unique(family_group()),
               col = unique(my_col),
               lty= 1,
               lwd = 5,
               cex=.7
        )
        dev.off()
        fs <- c(fs, paste("Heatmap_df1.", input$selectModule11_p5 ,"_df2.", input$selectModule22_p5,".pdf", sep = ""))
      }

    }else{
      svg(paste("Co_inertia_mainPlot_df1.", input$selectModule1_p5 ,"_df2.", input$selectModule2_p5,".svg", sep = ""), width = 10, height = 10)
      print(ggplot(data = selected_sampleInfo_all()) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], col = colnames(selected_sampleInfo_all())[3]), size = 3) +
              geom_point(aes_string(x = colnames(selected_sampleInfo_all())[4], y = colnames(selected_sampleInfo_all())[5], col = colnames(selected_sampleInfo_all())[3]), shape = 5, size = 3) +
              scale_shape(solid = TRUE) +
              geom_segment(aes_string(x = colnames(selected_sampleInfo_all())[1], y = colnames(selected_sampleInfo_all())[2], colour = colnames(selected_sampleInfo_all())[3], xend = colnames(selected_sampleInfo_all())[4], yend = colnames(selected_sampleInfo_all())[5])) +
              geom_label_repel(data = selected_axis1_drivers(),
                               aes(x = 0.10 * SV1, y = 0.10 * SV2, label = Spec, fill = feature_type_axis1),
                               size = 4, segment.size = 0.5,
                               label.padding = unit(0.1, "lines"), label.size = 0.5) +
              labs(x = sprintf("Axis1 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[1] / sum(mcoia()[["mcoa"]]$pseudoeig), 2)),
                   y = sprintf("Axis2 [%s%% Variance]",
                               100 * round(mcoia()[["mcoa"]]$pseudoeig[2] / sum(mcoia()[["mcoa"]]$pseudoeig), 2))))
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

      if ((input$selectModule11_p5 == "Individual_Variables" && is.null(input$selectVariablesExprDat) || (input$selectModule22_p5 == "Individual_Variables" && is.null(input$selectVariables)) )){
        print("something")

      }else{
        svg(paste("Heatmap_df1.", input$selectModule11_p5 ,"_df2.", input$selectModule22_p5,".svg", sep = ""), 10, 10)
        my_group=as.numeric(as.factor(family_group()))
        my_col=brewer.pal(9, "Set1")[my_group]
        heatmap(corr_expr(), ColSideColors = my_col, margins = c(15,15))

        legend("topright",
               legend = unique(family_group()),
               col = unique(my_col),
               lty= 1,
               lwd = 5,
               cex=.7
        )
        dev.off()
        fs <- c(fs, paste("Heatmap_df1.", input$selectModule11_p5 ,"_df2.", input$selectModule22_p5,".svg", sep = ""))
      }

    }

    zip(zipfile=filename, files=fs)


  },
  contentType = "application/zip")

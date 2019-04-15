################################################
#### SERVER FUNCTIONS : NETWORK EXPLORATION ####
################################################
# This script gathers all the functions related to the network exploration tabulation of MiBiOmics
# The first functions use the WGCNA functionality to show you the modules associations of your 
# network. The rest of the functions allow you to explore the specificity of your modules (the 
# relativa abundance if you are working with microbiote data, a representation of the module as a
# hive plot). A PLS can also be realized on each module to verify their association to an
# external trait and the content of the module can be downloaded.

#### REACTIVE OBJECTS ####

module_Relative_abundance <- reactive({
  modGenes = (selectedDynamicColor()== input$selectModule)
  exprDat <- exprDat_WGCNA()[which(rownames(exprDat_WGCNA()) %in% rownames(sampleAnnot_2())),]
  exprDat <- exprDat[match(rownames(exprDat), rownames(sampleAnnot_2())),]
  if (input$LoadExample == "No" && input$LoadExample2 == "No"){
    if (input$TypeAnalysis == "multivariate"){
      validate(
        #need(input$CountingT == "OTUs", "Only for OTUs counting table (multi-omics)."),
        need(input$TaxonFile1, "Only for OTUs counting table.(multi-omics)")
      )
    }else{
      validate(
        #need(input$CountingT1 == "OTUs", "Only for OTUs counting table (simple)."),
        need(input$TaxonFile, "The taxon annotation file was not uploaded.(simple)")
      )
    }
  }
  
  
  if (input$selectVariable_RelAb != "unordered"){
    exprDat <- exprDat[order(sampleAnnot_2()[,input$selectVariable_RelAb]),]
  }
  
  selectedExpr <- exprDat[, modGenes]
  if (input$Transformation =="Yes"){
    if (input$TypeTransformation=="CLR"){
      exprDat_rel_ab <- clrInv(selectedExpr)
    }else{
      exprDat_rel_ab <- as.data.frame(t(apply(selectedExpr, 1, function(x) x / sum(x))))
    }
  }else{
    exprDat_rel_ab <- as.data.frame(t(apply(selectedExpr, 1, function(x) x / sum(x))))
    
  }
  exprDat_rel_ab <- setDT(as.data.frame(t(exprDat_rel_ab)), keep.rownames = TRUE)
  exprDat_rel_ab <- merge(exprDat_rel_ab, taxTable1(), by = "rn")
  exprDat_rel_ab$rn <- as.character(exprDat_rel_ab$rn)
  exprDat_rel_ab$rn <- factor(exprDat_rel_ab$rn, levels = unique(exprDat_rel_ab$rn))
  df_long <- melt(exprDat_rel_ab)
  df_long
})

module_Relative_abundanceSec <- reactive({
  exprDat <- exprDatSec_WGCNA()
  
  modGenes = (selectedDynamicColor2() == input$selectModuleSec)
  
  
  if (input$selectVariable_RelAbSec != "unordered"){
    exprDat <- exprDat[order(sampleAnnot_sec()[,input$selectVariable_RelAbSec]),]
  }
  
  selectedExpr <- exprDat[, modGenes]
  if (input$Transformation1 =="Yes"){
    if (input$TypeTransformation1=="CLR"){
      exprDat_rel_ab <- clrInv(selectedExpr)
    }else{
      exprDat_rel_ab <- as.data.frame(t(apply(selectedExpr, 1, function(x) x / sum(x))))
    }
  }else{
    exprDat_rel_ab <- as.data.frame(t(apply(selectedExpr, 1, function(x) x / sum(x))))
  }
  exprDat_rel_ab <- setDT(as.data.frame(t(exprDat_rel_ab)), keep.rownames = TRUE)
  exprDat_rel_ab <- merge(exprDat_rel_ab, taxTable1(), by = "rn")
  exprDat_rel_ab$rn <- as.character(exprDat_rel_ab$rn)
  exprDat_rel_ab$rn <- factor(exprDat_rel_ab$rn, levels = unique(exprDat_rel_ab$rn))
  df_long <- melt(exprDat_rel_ab)
  df_long
})

# Order the MEs according to a variable inputed by the user 
MEs_ordered <- reactive ({
  MEs <- data.frame(selectedMEs(), rownames(sampleAnnot_2()))
  colnames(MEs) <- c(colnames(selectedMEs()), "sampleName")
  if (input$selectVariable_barplot != "unordered"){
    MEs <- MEs[order(sampleAnnot_2()[,input$selectVariable_barplot]),]
  }
  MEs$sampleName <- as.character(MEs$sampleName)
  MEs$sampleName <- factor(MEs$sampleName, levels = unique(MEs$sampleName))
  MEs
})

MEs_orderedSec <- reactive ({
  MEs <- data.frame(selectedMEs2(), rownames(sampleAnnot_sec()))
  colnames(MEs) <- c(colnames(selectedMEs2()), "sampleName")
  if (input$selectVariable_barplotSec != "unordered"){
    MEs <- MEs[order(sampleAnnot_sec()[,input$selectVariable_barplotSec]),]
  }
  MEs$sampleName <- as.character(MEs$sampleName)
  MEs$sampleName <- factor(MEs$sampleName, levels = unique(MEs$sampleName))
  MEs
})

sampleAnnot_ordered <- reactive({
  sampleAnnot <- data.frame(sampleAnnot_2(), rownames(sampleAnnot_2()))
  colnames(sampleAnnot) <- c(colnames(sampleAnnot_2()), "sampleName")
  if (input$selectVariable_barplot != "unordered"){
    sampleAnnot <- sampleAnnot[order(sampleAnnot_2()[,input$selectVariable_barplot]),]
  }
  sampleAnnot$sampleName <- as.character(sampleAnnot$sampleName)
  sampleAnnot$sampleName <- factor(sampleAnnot$sampleName, levels = unique(sampleAnnot$sampleName))
  sampleAnnot
})

sampleAnnot_orderedSec <- reactive({
  sampleAnnot <- data.frame(sampleAnnot_sec(), rownames(sampleAnnot_sec()))
  colnames(sampleAnnot) <- c(colnames(sampleAnnot_sec()), "sampleName")
  if (input$selectVariable_barplotSec != "unordered"){
    sampleAnnot <- sampleAnnot[order(sampleAnnot_sec()[,input$selectVariable_barplotSec]),]
  }
  sampleAnnot$sampleName <- as.character(sampleAnnot$sampleName)
  sampleAnnot$sampleName <- factor(sampleAnnot$sampleName, levels = unique(sampleAnnot$sampleName))
  sampleAnnot
})

geneModuleMembership <- reactive({
  modNames = substring(names(selectedMEs()), 3)
  geneModuleMembership = as.data.frame(cor(exprDat_WGCNA(), selectedMEs(), method= input$selectCorrelation, use="na.or.complete"));
  
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  geneModuleMembership
})

geneModuleMembershipSec <- reactive({
  modNames = substring(names(selectedMEs2()), 3)
  geneModuleMembership = as.data.frame(cor(exprDatSec_WGCNA(), selectedMEs2(), method= input$selectCorrelationSec, use="na.or.complete"));
  
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  geneModuleMembership
})

geneTraitSignificance <- reactive({
  colnames_used <- c()
  for (i in names(sampleAnnot_2())) {
    if(length(unique(sampleAnnot_2()[,i])) > 1){
      if (is.numeric(sampleAnnot_2()[,i])){
        env = as.data.frame(sampleAnnot_2()[,i]) 
      }else{
        env = as.data.frame(as.numeric(as.factor(sampleAnnot_2()[,i])))
      }
      colnames(env) = i
      colnames_used <- c(colnames_used, i)
      modNames = substring(names(selectedMEs()), 3)
      if (i == colnames(sampleAnnot_2())[1]){
        geneTraitSignificance = as.data.frame(cor(exprDat_WGCNA(), env, method= input$selectCorrelation, use="na.or.complete"));
      }else{
        geneTraitSignificance = data.frame(geneTraitSignificance, cor(exprDat_WGCNA(), env, method= input$selectCorrelation, use="na.or.complete"));
      }
    }
  }
  colnames(geneTraitSignificance) <- paste("GS.", colnames_used, sep = "")
  geneTraitSignificance
})

geneTraitSignificanceSec <- reactive({
  colnames_used <- c()
  for (i in names(sampleAnnot_sec())) {
    if(length(unique(sampleAnnot_sec()[,i])) > 1){
      if (is.numeric(sampleAnnot_sec()[,i])){
        env = as.data.frame(sampleAnnot_sec()[,i])
      }else{
        env = as.data.frame(as.numeric(as.factor(sampleAnnot_sec()[,i])))
      }
      colnames(env) = i
      colnames_used <- c(colnames_used, i)
      modNames = substring(names(selectedMEs2()), 3)
      if (i == colnames(sampleAnnot_sec())[1]){
        geneTraitSignificance = as.data.frame(cor(exprDatSec_WGCNA(), env, method= input$selectCorrelationSec, use="na.or.complete"));
      }else{
        geneTraitSignificance = data.frame(geneTraitSignificance, cor(exprDatSec_WGCNA(), env, method= input$selectCorrelationSec, use="na.or.complete"));
      }
    }
  }
  colnames(geneTraitSignificance) <- paste("GS.", colnames_used, sep = "")
  geneTraitSignificance
})


moduleMembership <- reactive ({
  module = input$selectModule
  column = paste("MM",module, sep = "")
  column2 = paste("GS.", input$selectVariable1, sep = "")
  moduleGenes = (selectedDynamicColor()==module)
  a <- abs(geneModuleMembership()[moduleGenes, column])
  b <- abs(geneTraitSignificance()[moduleGenes, column2])
  c <- data.frame(a,b)
  colnames(c) <- c("moduleMembership", "TraitSignificance")
  c
})

moduleMembershipSec <- reactive ({
  module = input$selectModuleSec
  column = paste("MM",module, sep = "")
  column2 = paste("GS.", input$selectVariable1Sec, sep = "")
  moduleGenes = (selectedDynamicColor2()==module)
  a <- abs(geneModuleMembershipSec()[moduleGenes, column])
  b <- abs(geneTraitSignificanceSec()[moduleGenes, column2])
  c <- data.frame(a,b)
  colnames(c) <- c("moduleMembership", "TraitSignificance")
  c
})

textModuleMembership <- reactive({
  module = input$selectModule
  column = paste("MM",module, sep = "")
  column2 = paste("GS.", input$selectVariable1, sep = "")
  moduleGenes = (selectedDynamicColor()== input$selectModule)
  pv = corPvalueStudent(cor(abs(geneModuleMembership()[moduleGenes, column]), abs(geneTraitSignificance()[moduleGenes, column2]), method= input$selectCorrelation, use="na.or.complete"), length(geneTraitSignificance()[moduleGenes, column2]))
  Correlation <- cor(abs(geneModuleMembership()[moduleGenes, column]), abs(geneTraitSignificance()[moduleGenes, column2]))
  text_scatter = paste("Corr = ", round(Correlation, 2), "\nP-val = ", pv)
  text_scatter
})


textModuleMembershipSec <- reactive({
  module = input$selectModuleSec
  column = paste("MM",module, sep = "")
  column2 = paste("GS.", input$selectVariable1Sec, sep = "")
  moduleGenes = (selectedDynamicColor2()== input$selectModuleSec)
  pv = corPvalueStudent(cor(abs(geneModuleMembershipSec()[moduleGenes, column]), abs(geneTraitSignificanceSec()[moduleGenes, column2]), method= input$selectCorrelationSec, use="na.or.complete"), length(geneTraitSignificanceSec()[moduleGenes, column2]))
  Correlation <- cor(abs(geneModuleMembershipSec()[moduleGenes, column]), abs(geneTraitSignificanceSec()[moduleGenes, column2]))
  text_scatter = paste("Corr = ", round(Correlation, 2), "\nP-val = ", pv)
  text_scatter
})


# Recalculate MEs with color labels
selectedCondDF <- reactive({
  moduleColors = selectedDynamicColor()
  nSamples <- nrow(exprDat_WGCNA())
  MEs0 = moduleEigengenes(exprDat_WGCNA(), moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  text_matrix_vector <- c()
  
  for (condition in 1:(ncol(sampleAnnot_2()))){
    
    condition_name <- colnames(sampleAnnot_2())[condition]
    if (is.numeric(sampleAnnot_2()[,condition_name])){
      moduleTraitCor = cor(selectedMEs(), sampleAnnot_2()[,condition_name], use = "p", method = input$selectCorrelation)
    }else{
      sampleAnnot4 <- sampleAnnot_2()
      sampleAnnot4[,condition_name] <- as.factor(sampleAnnot4[, condition_name])
      sampleAnnot4[,condition_name] <- as.numeric(sampleAnnot4[,condition_name])
      moduleTraitCor = cor(selectedMEs(), sampleAnnot4[,condition_name], use = "p", method = input$selectCorrelation)
    }
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
    colnames(moduleTraitCor) <- c(condition_name)
    # Will display correlations and their p-values
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "")
    text_matrix_vector <- c(text_matrix_vector, textMatrix)
    dim(textMatrix) = dim(moduleTraitCor)
    
    # create dataframe containing the correlation and the color
    if(condition == 1){
      
      Condition_df <- data.frame(moduleTraitCor)
      
    }else{
      
      Condition_df <- data.frame(Condition_df, moduleTraitCor)
      
    }
    
    melted_moduleTraitCor <- melt(moduleTraitCor)
    melted_moduleTraitCor <- data.frame(melted_moduleTraitCor, as.data.frame(substr(names(selectedMEs()), 3, 40)))
    colnames(melted_moduleTraitCor)[4] <- "color"
    
  }
  
  
  Condition_df <- setDT(Condition_df, keep.rownames = TRUE)[]
  melted_condition_df <- melt(Condition_df, id = 1)
  melted_condition_df$label <- as.data.frame(text_matrix_vector)
  melted_condition_df
  
})

# Recalculate MEs with color labels
selectedCondDFSec <- reactive({
  moduleColors = selectedDynamicColor2()
  nSamples <- nrow(exprDatSec_WGCNA())
  MEs0 = moduleEigengenes(exprDatSec_WGCNA(), moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  text_matrix_vector <- c()
  
  for (condition in 1:(ncol(sampleAnnot_sec()))){
    
    condition_name <- colnames(sampleAnnot_sec())[condition]
    if (is.numeric(sampleAnnot_sec()[,condition_name])){
      moduleTraitCor = cor(selectedMEs2(), sampleAnnot_sec()[,condition_name], use = "p", method = input$selectCorrelationSec)
    }else{
      sampleAnnot4 <- sampleAnnot_sec()
      sampleAnnot4[,condition_name] <- as.factor(sampleAnnot4[, condition_name])
      sampleAnnot4[,condition_name] <- as.numeric(sampleAnnot4[,condition_name])
      moduleTraitCor = cor(selectedMEs2(), sampleAnnot4[,condition_name], use = "p", method = input$selectCorrelationSec)
    }
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
    colnames(moduleTraitCor) <- c(condition_name)
    # Will display correlations and their p-values
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "")
    text_matrix_vector <- c(text_matrix_vector, textMatrix)
    dim(textMatrix) = dim(moduleTraitCor)
    
    # create dataframe containing the correlation and the color
    if(condition == 1){
      
      Condition_df <- data.frame(moduleTraitCor)
      
    }else{
      
      Condition_df <- data.frame(Condition_df, moduleTraitCor)
      
    }
    
    melted_moduleTraitCor <- melt(moduleTraitCor)
    melted_moduleTraitCor <- data.frame(melted_moduleTraitCor, as.data.frame(substr(names(selectedMEs2()), 3, 40)))
    colnames(melted_moduleTraitCor)[4] <- "color"
    
  }
  
  Condition_df <- setDT(Condition_df, keep.rownames = TRUE)[]
  melted_condition_df <- melt(Condition_df, id = 1)
  melted_condition_df$label <- as.data.frame(text_matrix_vector)
  melted_condition_df
  
})


#### VIP SCORE

markers.data <- reactive({
  markers <- colnames(exprDat_WGCNA())[selectedDynamicColor()==input$colorModule]
  exprDat_WGCNA()[,markers]
  
})

env.markers.data <- reactive({
  sampleAnnot <- sampleAnnot_2()
  if (!is.numeric(sampleAnnot[,input$sampleAnnotSelection])){
    sampleAnnot[,input$sampleAnnotSelection] <- as.numeric(as.factor(sampleAnnot[,input$sampleAnnotSelection]))
  }
  env = as.data.frame(sampleAnnot[,input$sampleAnnotSelection])
  colnames(env) = input$sampleAnnotSelection
  data.frame(cbind(markers.data(), env))
})

fit.object.ncomp <- reactive({
  
  if (length(colnames(exprDat_WGCNA())[selectedDynamicColor()==input$colorModule]) >= 10){
    ncomp=10
  }else{
    ncomp=length(colnames(exprDat_WGCNA())[selectedDynamicColor()==input$colorModule])
  }
  formulaChar <- paste(input$sampleAnnotSelection, "~.", sep ="")
  fit = plsr(formula(formulaChar), data=env.markers.data(), validation="LOO", method = "oscorespls", ncomp=ncomp)
  
  fit
  
})

fit.object <- reactive({
  ncomp=input$ncomponent
  formulaChar <- paste(input$sampleAnnotSelection, "~.", sep ="")
  fit = plsr(formula(formulaChar), data=env.markers.data(), validation="LOO", method = "oscorespls", ncomp=ncomp)
  fit
})


cvs.object <- reactive({
  ncomp = input$ncomponent
  fit <- fit.object()
  env.markers.data <- env.markers.data()
  cvs = round( cor(fit$val$pred[,,ncomp], env.markers.data[[input$sampleAnnotSelection]][!is.na(env.markers.data[[input$sampleAnnotSelection]])]), digits=3);
  cvs
})

vip.df <- reactive({
  
  vip = vector("numeric", ncol(markers.data())); for(j in 1:ncol(markers.data()))  vip[j] <- VIPjh(fit.object(), j, input$ncomponent)
  vip = round(vip, digits=2)
  o = order(vip, decreasing=TRUE)
  df <- data.frame(Marker = colnames(markers.data()), VIP = vip, Temp.Cor = geneTraitSignificance()[colnames(markers.data()),])
  df
})

#### HIVE PLOT

edge_to_node <- reactive({
  module <- input$colorModule
  probes <- colnames(exprDat_WGCNA())
  inModule <- is.finite(match(selectedDynamicColor(), module))
  modProbes <- probes[inModule]
  modTOM <- selectedDissTOM()[inModule, inModule]
  dimnames(modTOM) <- list(modProbes, modProbes)
  cyt <- exportNetworkToCytoscape(modTOM,
                                  edgeFile = NULL,
                                  nodeFile = NULL,
                                  weighted = TRUE,
                                  nodeNames = modProbes,
                                  nodeAttr = selectedDynamicColor()[inModule])
  cyt
})

df.hive.plot <- reactive({
  vip.df <- vip.df()
  col <- paste("Temp.Cor.GS.", input$sampleAnnotSelection, sep = "")
  if (input$TaxonFile || input$TaxonFile1 || input$LoadExample == "Yes" || input$LoadExample2 == "Yes"){
    Annot <- taxTable1()[rownames(vip.df), input$taxAnnot]
  }else{
    Annot <- rownames(vip.df)
  }
  df.hive.plot <- data.frame(label = rownames(vip.df), x1 = vip.df[["VIP"]], x2 = rep(0, nrow(vip.df)), y1 = rep(0, nrow(vip.df)), y2 = vip.df[[col]], annotation = Annot)
  df.hive.plot
})

df.hive.to.plot <- reactive({
  edge_node <- edge_to_node()
  edgeData <- edge_node$`edgeData`
  vector_geom_curve(edgeData, df.hive.plot())
})

#### VIP SCORE

markers.data_D2 <- reactive({
  markers <- colnames(exprDatSec_WGCNA())[selectedDynamicColor2()==input$colorModule_D2]
  exprDatSec_WGCNA()[,markers]
})

env.markers.data_D2 <- reactive({
  sampleAnnot <- sampleAnnot_sec()
  if (!is.numeric(sampleAnnot[,input$sampleAnnotSelection_D2])){
    sampleAnnot[,input$sampleAnnotSelection_D2] <- as.numeric(as.factor(sampleAnnot[,input$sampleAnnotSelection_D2]))
  }
  env = as.data.frame(sampleAnnot[,input$sampleAnnotSelection_D2])
  colnames(env) = input$sampleAnnotSelection_D2
  data.frame(cbind(markers.data_D2(), env))
})

fit.object.ncomp_D2 <- reactive({
  if (length(colnames(exprDatSec_WGCNA())[selectedDynamicColor2()==input$colorModule_D2]) >= 10){
    ncomp=10
  }else{
    ncomp=length(colnames(exprDatSec_WGCNA())[selectedDynamicColor2()==input$colorModule_D2])
  }
  formulaChar <- paste(input$sampleAnnotSelection_D2, "~.", sep ="")
  fit = plsr(formula(formulaChar), data=env.markers.data_D2(), validation="LOO", method = "oscorespls", ncomp=ncomp)
  fit
})

fit.object_D2 <- reactive({
  ncomp=input$ncomponent_D2
  formulaChar <- paste(input$sampleAnnotSelection_D2, "~.", sep ="")
  fit = plsr(formula(formulaChar), data=env.markers.data_D2(), validation="LOO", method = "oscorespls", ncomp=ncomp)
  fit
})


cvs.object_D2 <- reactive({
  ncomp = input$ncomponent_D2
  fit <- fit.object_D2()
  env.markers.data <- env.markers.data_D2()
  cvs = round( cor(fit$val$pred[,,ncomp], env.markers.data[[input$sampleAnnotSelection_D2]][!is.na(env.markers.data[[input$sampleAnnotSelection_D2]])]), digits=3);
  cvs
})

vip.df_D2 <- reactive({
  
  vip = vector("numeric", ncol(markers.data_D2())); for(j in 1:ncol(markers.data_D2()))  vip[j] <- VIPjh(fit.object_D2(), j, input$ncomponent_D2)
  vip = round(vip, digits=2)
  o = order(vip, decreasing=TRUE)
  df <- data.frame(Marker = colnames(markers.data_D2()), VIP = vip, Temp.Cor = geneTraitSignificanceSec()[colnames(markers.data_D2()),])
  df
})

#### HIVE PLOT

edge_to_node_D2 <- reactive({
  module <- input$colorModule_D2
  probes <- colnames(exprDatSec_WGCNA())
  inModule <- is.finite(match(selectedDynamicColor2(), module))
  modProbes <- probes[inModule]
  modTOM <- selectedDissTOM2()[inModule, inModule]
  dimnames(modTOM) <- list(modProbes, modProbes)
  cyt <- exportNetworkToCytoscape(modTOM,
                                  edgeFile = NULL,
                                  nodeFile = NULL,
                                  weighted = TRUE,
                                  nodeNames = modProbes,
                                  nodeAttr = selectedDynamicColor2()[inModule])
  cyt
})

df.hive.plot_D2 <- reactive({
  vip.df <- vip.df_D2()
  col <- paste("Temp.Cor.GS.", input$sampleAnnotSelection_D2, sep = "")
  if (input$LoadExample2 == "Yes"){
    Annot <- rownames(vip.df)
    
  }else{
    if (input$OmicTable == "OTUs" && input$TaxonFile1){
      Annot <- taxTable1()[rownames(vip.df), input$taxAnnotSec]
    }else{
      Annot <- rownames(vip.df)
    }
  }
  df.hive.plot <- data.frame(label = rownames(vip.df), x1 = vip.df[["VIP"]], x2 = rep(0, nrow(vip.df)), y1 = rep(0, nrow(vip.df)), y2 = vip.df[[col]], annotation = Annot)
  df.hive.plot
})

df.hive.to.plot_D2 <- reactive({
  edge_node <- edge_to_node_D2()
  edgeData <- edge_node$`edgeData`
  vector_geom_curve(edgeData, df.hive.plot_D2())
})

#### INTERFACE VARIABLES ####

output$SelectVariable1 <- renderUI({
  selectInput("selectVariable1", 
              label = "Choose a variable: ", 
              choices = colnames(sampleAnnot_2()), 
              selected = colnames(sampleAnnot_2())[2])
})

output$SelectVariable1Sec <- renderUI({
  selectInput("selectVariable1Sec", 
              label = "Choose a variable: ", 
              choices = colnames(sampleAnnot_2()), 
              selected = colnames(sampleAnnot_2())[2])
})

output$sampleAnnotSelection <- renderUI({
  selectInput("sampleAnnotSelection",
              label = "Choose a variable: ",
              choices = colnames(sampleAnnot_2()),
              selected = colnames(sampleAnnot_2())[2])
})

output$colorModule <-  renderUI({
  selectInput("colorModule",
              label= "Choose a module color: ",
              choices = substr(names(selectedMEs()), 3, 40))
})


output$sampleAnnotSelection_D2 <- renderUI({
  selectInput("sampleAnnotSelection_D2",
              label = "Choose a variable: ",
              choices = colnames(sampleAnnot_2()),
              selected = colnames(sampleAnnot_2())[2])
})

output$colorModule_D2 <-  renderUI({
  selectInput("colorModule_D2",
              label= "Choose a module color: ",
              choices = substr(names(selectedMEs2()), 3, 40))
})

output$SelectVariable_barplot <- renderUI({
  selectInput("selectVariable_barplot", 
              label = "Choose a variable: ", 
              choices = c("unordered", colnames(sampleAnnot_2())), 
              selected = "unordered")
})

output$SelectVariable_barplotSec <- renderUI({
  selectInput("selectVariable_barplotSec", 
              label = "Choose a variable: ", 
              choices = c("unordered", colnames(sampleAnnot_2())), 
              selected = "unordered")
})

output$SelectVariable_RelAb <- renderUI({
  selectInput("selectVariable_RelAb", 
              label = "Choose a variable: ", 
              choices = c("unordered", colnames(sampleAnnot_2())),  
              selected = "unordered")
})

output$SelectVariable_RelAbSec <- renderUI({
  selectInput("selectVariable_RelAbSec", 
              label = "Choose a variable: ", 
              choices = c("unordered", colnames(sampleAnnot_2())),  
              selected = "unordered")
})

output$SelectModule <- renderUI({
  selectInput("selectModule",
              label= "Choose a module color: ",
              choices = substr(names(selectedMEs()), 3, 40))
})

output$SelectModuleSec <- renderUI({
  selectInput("selectModuleSec",
              label= "Choose a module color: ",
              choices = substr(names(selectedMEs2()), 3, 40))
})

output$SelectTaxo1 <- renderUI({
  selectInput("selectTaxo2",
              label = "Select a taxonomic level for the relative abundance plot",
              choices = colnames(taxTable_present()),
              selected = colnames(taxTable_present())[3])
})

output$taxAnnot <- renderUI({
  selectInput("taxAnnot",
              label = "Select a taxonomic level for the relative abundance plot",
              choices = colnames(taxTable_present()),
              selected = colnames(taxTable_present())[3])
})

output$taxAnnotSec <- renderUI({
  selectInput("taxAnnotSec",
              label = "Select a taxonomic level for the relative abundance plot",
              choices = colnames(taxTable_present()),
              selected = colnames(taxTable_present())[3])
})

output$SelectTaxo1Sec <- renderUI({
  selectInput("selectTaxo2Sec",
              label = "Select a taxonomic level for the relative abundance plot",
              choices = colnames(taxTable_present()),
              selected = colnames(taxTable_present())[3])
})

#### PLOT OUTPUTS ####


# Contribution of the samples to each module
output$Sample_Contribution <- renderPlot({
  moduleName <- input$selectModule
  MODS <- paste("ME", input$selectModule, sep = "")
  
  MEs_barplot <- 
    ggplot(data = MEs_ordered(), aes_string(x = "sampleName", y = MODS )) +
    geom_bar(stat = "identity", fill = moduleName) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  if (input$selectVariable_barplot != "unordered"){
    sampleCurve <- 
      ggplot(data = sampleAnnot_ordered(), aes_string(x = "sampleName", y = input$selectVariable_barplot, group = 1)) + 
      geom_point() + 
      theme_bw() +
      geom_line() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    plot_grid(MEs_barplot, sampleCurve, labels = c("A", "B"), nrow = 2, align = "v")
  }else{
    MEs_barplot
  }
})

# Contribution of the samples to each module dataset 2
output$Sample_ContributionSec <- renderPlot({
  moduleName <- input$selectModuleSec
  MODS <- paste("ME", input$selectModuleSec, sep = "")
  
  MEs_barplotSec <-
    ggplot(data = MEs_orderedSec(), aes_string(x = "sampleName", y = MODS )) +
    geom_bar(stat = "identity", fill = moduleName) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  if (input$selectVariable_barplotSec != "unordered"){
    sampleCurveSec <-
      ggplot(data = sampleAnnot_orderedSec(), aes_string(x = "sampleName", y = input$selectVariable_barplotSec, group = 1)) +
      geom_point() +
      theme_bw() +
      geom_line() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    plot_grid(MEs_barplotSec, sampleCurveSec, labels = c("A", "B"), nrow = 2, align = "v")
  }else{
    MEs_barplotSec
  }
})

# Relative abundance of the modules

output$Relative_Abundance_Module <- renderPlot({
  validate(
    need(!is.na(unique(module_Relative_abundance()[[input$selectTaxo2]])), "Uncharacterized at this taxonomic level")
  )
  ggplot(module_Relative_abundance(), aes_string(x = "variable", y = "value", fill = input$selectTaxo2)) +
    geom_bar(stat = "identity") +
    labs(x = "Samples", y = "Relative Abundance", title = paste("Relative abundance at the ", input$selectTaxo2, " Taxonomic level of the ", input$selectModule, " module.", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
})

output$Relative_Abundance_ModuleSec <- renderPlot({
  ggplot(module_Relative_abundanceSec(), aes_string(x = "variable", y = "value", fill = input$selectTaxo2Sec)) +
    geom_bar(stat = "identity") +
    labs(x = "Samples", y = "Relative Abundance", title = paste("Relative abundance at the ", input$selectTaxo2Sec, " Taxonomic level of the ", input$selectModuleSec, " module.", sep = "")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
})

# Relate Modules to External Trait

output$Corr_External_Trait <- renderPlot({
  ggplot(data = selectedCondDF(), aes(x = rn, y = variable, fill = value)) + 
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name=paste("Module\n Conditions \nCorrelation"))+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 10, hjust = 1))+
    geom_text(aes(label = label), color = "black", size = 3) +
    coord_fixed() +
    xlab('Module') +
    ylab('')
})


output$Corr_External_TraitSec <- renderPlot({
  ggplot(data = selectedCondDFSec(), aes(x = rn, y = variable, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name=paste("Module\n Conditions \nCorrelation"))+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 10, hjust = 1))+
    geom_text(aes(label = label), color = "black", size = 3) +
    coord_fixed() +
    xlab('Module') +
    ylab('')
})

# Module Membership Versus Correlation to External Trait

output$Module_Membership <- renderPlot({
  validate(
    need(length(unique(sampleAnnot_2()[, input$selectVariable1])) >1, "The external traits does not discriminate several conditions" )
  )
  ggplot(moduleMembership(), aes(x = moduleMembership, y = TraitSignificance)) +
    geom_point(color = input$selectModule) +
    theme_bw() +
    geom_smooth(method=lm, color = "white") +
    labs(x = paste("Module Membership in", input$selectModule, "module"), y = paste("OTU significance for", input$selectVariable1),
         title = paste("Module membership vs. gene significance\n", textModuleMembership() ,sep = "")) +
    geom_density_2d(color = '#E69F00')
})

output$Module_MembershipSec <- renderPlot({
  validate(
    need(length(unique(sampleAnnot_2()[, input$selectVariable1])) >1, "The external traits does not discriminate several conditions" )
  )
  ggplot(moduleMembershipSec(), aes(x = moduleMembership, y = TraitSignificance)) +
    geom_point(color = input$selectModuleSec) +
    theme_bw() +
    geom_smooth(method=lm, color = "white") +
    labs(x = paste("Module Membership in", input$selectModuleSec, "module"), y = paste("OTU significance for", input$selectVariable1Sec),
         title = paste("Module membership vs. gene significance\n", textModuleMembershipSec() ,sep = "")) +
    geom_density_2d(color = '#E69F00')
})

#### VIP and PLS 

output$PLS <- renderPlot({
  predplot(fit.object(), main=paste("Cor. in LOO CV:",cvs.object()), col="blue", cex=0.5)
  plot(fit.object(), ncomp=input$ncomponent, asp=1, line=TRUE)
})

output$CVS <- renderText({
  print(cvs.object())
})

output$ncomp <- renderPlot({
  plot(RMSEP(fit.object.ncomp(), legendpos = "topright"))
})

output$summaryFit <- renderText({
  print(summary(fit.object()))
})

#### HIVE PLOT 
output$edge_node <- renderPlot({
  df_hive<-df.hive.to.plot()

  #edge_to_node()
  ggplot() +
    geom_point(data = df.hive.plot(), aes(x1, y1, color = annotation, size = 3)) + # x1 is the VIP
    geom_text(data = df.hive.plot(), vjust = 0.5, angle = 45, aes(x1, y1, label = annotation)) +
    geom_point(data = df.hive.plot(), aes(x2, abs(y2), color = annotation, size = 3)) + # y2 is the correlation to the pH at time of filtering
    geom_text(data = df.hive.plot(), hjust = 0.5, angle = 45, aes(x2,abs(y2), label = annotation)) +      
    geom_curve(data = df.hive.to.plot(), curvature = 0.2, color = "grey ", size = 0.2,aes(x = vector_x, y = vector_y, xend = vector_xend, yend = abs(vector_yend))) +
    theme_gdocs() +
    scale_colour_tableau(palette = "Tableau 20") +
    labs(x = "VIP", y = paste("Spearman Correlation to", input$sampleAnnotSelection), colour = "Annotation")
  
})


#### VIP and PLS 

output$PLS_D2 <- renderPlot({
  predplot(fit.object_D2(), main=paste("Cor. in LOO CV:",cvs.object_D2()), col="blue", cex=0.5)
  plot(fit.object_D2(), ncomp=input$ncomponent_D2, asp=1, line=TRUE)
})

output$CVS_D2 <- renderText({
  print(cvs.object_D2())
})

output$ncomp_D2 <- renderPlot({
  plot(RMSEP(fit.object.ncomp_D2(), legendpos = "topright"))
})

output$summaryFit_D2 <- renderText({
  print(summary(fit.object_D2()))
})

#### HIVE PLOT 
output$edge_node_D2 <- renderPlot({
  ggplot() +
    geom_point(data = df.hive.plot_D2(), aes(x1, y1, color = annotation, size = 3)) + # x1 is the VIP
    geom_text(data = df.hive.plot_D2(), vjust = 0.5, angle = 45, aes(x1, y1, label = annotation)) +
    geom_point(data = df.hive.plot_D2(), aes(x2, abs(y2), color = annotation, size = 3)) + # y2 is the correlation to the pH at time of filtering
    geom_text(data = df.hive.plot_D2(), hjust = 0.5, angle = 45, aes(x2, abs(y2), label = annotation)) + # y2 is the correlation to the pH at time of filtering
    geom_curve(data = df.hive.to.plot_D2(), curvature = 0.2, color = "grey ", size = 0.2,aes(x = vector_x, y = vector_y, xend = vector_xend, yend = abs(vector_yend))) +
    theme_gdocs() +
    scale_colour_tableau(palette = "Tableau 20") +
    labs(x = "VIP", y = paste("Spearman Correlation to", input$sampleAnnotSelection_D2), colour = "Annotations")
  
})



#### DOWNLOADS ####

output$Download_Network_Exploration <- downloadHandler(
  filename = function(){
    paste("NETWORK_EXPLORATION.zip")
  },
  content = function (filename){
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    if (input$pdf_or_svg_p4 == "pdf"){
      pdf("Contribution_MEs.pdf", width = input$widthPDF, height = input$heightPDF)
      for (i in 1:ncol(selectedMEs())){
        print(barplot(selectedMEs()[,i], col= substr(colnames(selectedMEs())[i], 3, 40), main="", cex.main=2,
                      ylab="eigengene expression", cex.names = 0.7, names.arg = rownames(sampleAnnot_2()), las = 2))
      }
      dev.off()
      fs <- c(fs, "Contribution_MEs.pdf")
      pdf(paste("Module_", input$selectModule ,"_Corr_Vs_Membership.pdf", sep=""), width = input$widthPDF, height = input$heightPDF)
      print(ggplot(moduleMembership(), aes(x = moduleMembership, y = TraitSignificance)) +
              geom_point(color = input$selectModule) +
              geom_smooth(method=lm, color = "white") +
              labs(x = paste("Module Membership in", input$selectModule, "module"), y = paste("OTU significance for", input$selectVariable1),
                   title = paste("Module membership vs. gene significance\n", textModuleMembership() ,sep = "")) +
              geom_density_2d(color = '#E69F00'))
      
      dev.off()
      fs <- c(fs, paste("Module_", input$selectModule ,"_Corr_Vs_Membership.pdf", sep=""))
      if (input$TaxonFile || input$LoadExample == "Yes" || input$LoadExample2 == "Yes"){
        pdf(paste("Module_", input$selectModule ,"_Relative_Abundance.pdf", sep=""), width = input$widthPDF, height = input$heightPDF)
        
        print(ggplot(module_Relative_abundance(), aes_string(x = "variable", y = "value", fill = input$selectTaxo2)) +
                geom_bar(stat = "identity") +
                labs(x = "Samples", y = "Relative Abundance", title = paste("Relative abundance at the ", input$selectTaxo2, " Taxonomic level of the ", input$selectModule, " module.", sep = "")) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)))
        
        dev.off()
        fs <- c(fs, paste("Module_", input$selectModule ,"_Relative_Abundance.pdf", sep=""))
      }
      
      pdf("Module_Correlation_External_Trait.pdf", width = input$widthPDF, height = input$heightPDF)
      print(ggplot(data = selectedCondDF(), aes(x = rn, y = variable, fill = value)) +
              geom_tile(color = "white") +
              scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                   midpoint = 0, limit = c(-1,1), space = "Lab",
                                   name=paste("Module\n Conditions \nCorrelation"))+
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                               size = 10, hjust = 1))+
              geom_text(aes(label = label), color = "black", size = 3) +
              coord_fixed() +
              xlab('Module') +
              ylab(''))
      dev.off()
      fs <- c(fs, "Module_Correlation_External_Trait.pdf")
    }else{
      for (i in 1:ncol(selectedMEs())){
        svg(paste("Contribution_MEs_", substr(colnames(selectedMEs())[i], 3, 40), ".svg", sep = ""), width = input$widthPDF, height = input$heightPDF)
        print(barplot(selectedMEs()[,i], col= substr(colnames(selectedMEs())[i], 3, 40), main="", cex.main=2,
                      ylab="eigengene expression", cex.names = 0.7, names.arg = rownames(sampleAnnot_2()), las = 2))
        dev.off()
        fs <- c(fs, paste("Contribution_MEs_", substr(colnames(selectedMEs())[i], 3, 40), ".svg", sep = ""))
      }
      svg(paste("Module_", input$selectModule ,"_Corr_Vs_Membership.svg", sep=""), width = input$widthPDF, height = input$heightPDF)
      print(ggplot(moduleMembership(), aes(x = moduleMembership, y = TraitSignificance)) +
              geom_point(color = input$selectModule) +
              geom_smooth(method=lm, color = "white") +
              labs(x = paste("Module Membership in", input$selectModule, "module"), y = paste("OTU significance for", input$selectVariable1),
                   title = paste("Module membership vs. gene significance\n", textModuleMembership() ,sep = "")) +
              geom_density_2d(color = '#E69F00'))
      
      dev.off()
      fs <- c(fs, paste("Module_", input$selectModule ,"_Corr_Vs_Membership.svg", sep=""))
      if (input$TaxonFile || input$LoadExample == "Yes" || input$LoadExample2 == "Yes"){
        svg(paste("Module_", input$selectModule ,"_Relative_Abundance.svg", sep=""), width = input$widthPDF, height = input$heightPDF)
        
        print(ggplot(module_Relative_abundance(), aes_string(x = "variable", y = "value", fill = input$selectTaxo2)) +
                geom_bar(stat = "identity") +
                labs(x = "Samples", y = "Relative Abundance", title = paste("Relative abundance at the ", input$selectTaxo2, " Taxonomic level of the ", input$selectModule, " module.", sep = "")) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)))
        
        dev.off()
        fs <- c(fs, paste("Module_", input$selectModule ,"_Relative_Abundance.svg", sep=""))
      }
      
      svg("Module_Correlation_External_Trait.svg", width = input$widthPDF, height = input$heightPDF)
      print(ggplot(data = selectedCondDF(), aes(x = rn, y = variable, fill = value)) +
              geom_tile(color = "white") +
              scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                   midpoint = 0, limit = c(-1,1), space = "Lab",
                                   name=paste("Module\n Conditions \nCorrelation"))+
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                               size = 10, hjust = 1))+
              geom_text(aes(label = label), color = "black", size = 3) +
              coord_fixed() +
              xlab('Module') +
              ylab(''))
      dev.off()
      fs <- c(fs, "Module_Correlation_External_Trait.svg")
    }
    
    zip(zipfile=filename, files=fs)
    
  },
  contentType = "application/zip")

#### PLS
output$PLS_VIP <- downloadHandler(
  filename = function(){
    paste(input$sampleAnnotSelection, "_", input$colorModule, "markers.vip.csv", sep="")
  },
  content = function (filename){
    vip <- vip.df()
    if (input$LoadExample == "Yes" || input$LoadExample2 == "Yes" || input$TaxonFile || input$TaxonFile1){
      taxTable <- taxTable1()[which(rownames(taxTable1()) %in% rownames(vip)),]
      taxTable <- taxTable[match(rownames(taxTable),rownames(vip)),]
      
      vip_taxAnnot <- merge(vip, taxTable, by=0, all=TRUE)
    }else{
      vip_taxAnnot <- vip
    }

    write.csv(vip_taxAnnot, filename)
  }
)

output$hivePlot_D1 <- downloadHandler(
  filename = function(){
    paste("HIVE.zip")
  },
  content = function (filename){
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    if (input$pdf_or_svg_hive_p4 == "pdf"){
      pdf("HIVE.pdf", width = input$HIVEwidthPDF, height = input$HIVEheightPDF)

        print(  ggplot() +
                  geom_point(data = df.hive.plot(), aes(x1, y1, color = annotation, size = 3)) + # x1 is the VIP
                  geom_text(data = df.hive.plot(), vjust = 0.5, angle = 45, aes(x1, y1, label = annotation)) +
                  geom_point(data = df.hive.plot(), aes(x2, abs(y2), color = annotation, size = 3)) + # y2 is the correlation to the pH at time of filtering
                  geom_text(data = df.hive.plot(), hjust = 0.5, angle = 45, aes(x2,abs(y2), label = annotation)) +      
                  geom_curve(data = df.hive.to.plot(), curvature = 0.2, color = "grey ", size = 0.2,aes(x = vector_x, y = vector_y, xend = vector_xend, yend = abs(vector_yend))) +
                  theme_gdocs() +
                  scale_colour_tableau(palette = "Tableau 20") +
                  labs(x = "VIP", y = paste("Spearman Correlation to", input$sampleAnnotSelection), colour = "Annotation"))
        dev.off()
        fs <- c(fs, "HIVE.pdf")
      
    }else{

        svg("HIVE.svg", width = input$HIVEwidthPDF, height = input$HIVEwidthPDF)
      
      print(  ggplot() +
                geom_point(data = df.hive.plot(), aes(x1, y1, color = annotation, size = 3)) + # x1 is the VIP
                geom_text(data = df.hive.plot(), vjust = 0.5, angle = 45, aes(x1, y1, label = annotation)) +
                geom_point(data = df.hive.plot(), aes(x2, abs(y2), color = annotation, size = 3)) + # y2 is the correlation to the pH at time of filtering
                geom_text(data = df.hive.plot(), hjust = 0.5, angle = 45, aes(x2,abs(y2), label = annotation)) +      
                geom_curve(data = df.hive.to.plot(), curvature = 0.2, color = "grey ", size = 0.2,aes(x = vector_x, y = vector_y, xend = vector_xend, yend = abs(vector_yend))) +
                theme_gdocs() +
                scale_colour_tableau(palette = "Tableau 20") +
                labs(x = "VIP", y = paste("Spearman Correlation to", input$sampleAnnotSelection), colour = "Annotation"))
      dev.off()
      fs <- c(fs, "HIVE.pdf")

      
    }
    
    zip(zipfile=filename, files=fs)
    
  },
  contentType = "application/zip")

#### PLS 
output$PLS_VIP_D2 <- downloadHandler(
  filename = function(){
    paste(input$sampleAnnotSelection_D2, "_", input$colorModule_D2, "markers.vip_D2.csv", sep="")
  },
  content = function (filename){
    
    vip <- vip.df_D2()
    if (input$TaxonFile1){
      taxTable <- taxTable1()[which(rownames(taxTable1()) %in% rownames(vip)),]
      taxTable <- taxTable[match(rownames(taxTable),rownames(vip)),]
      
      vip_taxAnnot <- merge(vip, taxTable, by=0, all=TRUE)
    }else{
      vip_taxAnnot <- vip
    }
    
    write.csv(vip_taxAnnot, filename)
  }
)

output$hivePlot_D2 <- downloadHandler(
  filename = function(){
    paste("HIVE_D2.zip")
  },
  content = function (filename){
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    if (input$pdf_or_svg_hiveD2_p4 == "pdf"){
      pdf("HIVE_D2.pdf", width = input$HIVED2widthPDF, height = input$HIVED2heightPDF)
      
      print(  ggplot() +
                geom_point(data = df.hive.plot(), aes(x1, y1, color = annotation, size = 3)) + # x1 is the VIP
                geom_text(data = df.hive.plot(), vjust = 0.5, angle = 45, aes(x1, y1, label = annotation)) +
                geom_point(data = df.hive.plot(), aes(x2, abs(y2), color = annotation, size = 3)) + # y2 is the correlation to the pH at time of filtering
                geom_text(data = df.hive.plot(), hjust = 0.5, angle = 45, aes(x2,abs(y2), label = annotation)) +      
                geom_curve(data = df.hive.to.plot(), curvature = 0.2, color = "grey ", size = 0.2,aes(x = vector_x, y = vector_y, xend = vector_xend, yend = abs(vector_yend))) +
                theme_gdocs() +
                scale_colour_tableau(palette = "Tableau 20") +
                labs(x = "VIP", y = paste("Spearman Correlation to", input$sampleAnnotSelection), colour = "Annotation"))
      dev.off()
      fs <- c(fs, "HIVE_D2.pdf")
      
    }else{
      
      svg("HIVE_D2.svg", width = input$HIVED2widthPDF, height = input$HIVED2widthPDF)
      
      print(  ggplot() +
                geom_point(data = df.hive.plot(), aes(x1, y1, color = annotation, size = 3)) + # x1 is the VIP
                geom_text(data = df.hive.plot(), vjust = 0.5, angle = 45, aes(x1, y1, label = annotation)) +
                geom_point(data = df.hive.plot(), aes(x2, abs(y2), color = annotation, size = 3)) + # y2 is the correlation to the pH at time of filtering
                geom_text(data = df.hive.plot(), hjust = 0.5, angle = 45, aes(x2,abs(y2), label = annotation)) +      
                geom_curve(data = df.hive.to.plot(), curvature = 0.2, color = "grey ", size = 0.2,aes(x = vector_x, y = vector_y, xend = vector_xend, yend = abs(vector_yend))) +
                theme_gdocs() +
                scale_colour_tableau(palette = "Tableau 20") +
                labs(x = "VIP", y = paste("Spearman Correlation to", input$sampleAnnotSelection), colour = "Annotation"))
      dev.off()
      fs <- c(fs, "HIVE_D2.pdf")
      
      
    }
    
    zip(zipfile=filename, files=fs)
    
  },
  contentType = "application/zip")



#### DataSet 2 

output$Download_Network_Exploration_dataset_2 <- downloadHandler(
  filename = function(){
    paste("NETWORK_EXPLORATION_dataset2.zip")
  },
  content = function (filename){
    fs <- c()
    tmpdir <- tempdir()
    setwd(tempdir())
    if (input$pdf_or_svg_p4_dataset2 == "pdf"){
      pdf("Contribution_MEs.pdf", width = input$widthPDF_D2, height = input$heightPDF_D2)
      for (i in 1:ncol(selectedMEs2())){
        print(barplot(selectedMEs2()[,i], col= substr(colnames(selectedMEs2())[i], 3, 40), main="", cex.main=2,
                      ylab="eigengene expression", cex.names = 0.7, names.arg = rownames(sampleAnnot_sec()), las = 2))
      }
      dev.off()
      fs <- c(fs, "Contribution_MEs.pdf")
      
      pdf(paste("Module_", input$selectModuleSec ,"_Corr_Vs_Membership.pdf", sep=""), width = input$widthPDF_D2, height = input$heightPDF_D2)
      print(ggplot(moduleMembershipSec(), aes(x = moduleMembership, y = TraitSignificance)) +
              geom_point(color = input$selectModuleSec) +
              geom_smooth(method=lm, color = "white") +
              labs(x = paste("Module Membership in", input$selectModuleSec, "module"), y = paste("OTU significance for", input$selectVariable1Sec),
                   title = paste("Module membership vs. gene significance\n", textModuleMembershipSec() ,sep = "")) +
              geom_density_2d(color = '#E69F00'))
      
      dev.off()
      fs <- c(fs, paste("Module_", input$selectModuleSec ,"_Corr_Vs_Membership.pdf", sep=""))
      if (input$TaxonFile1){
        pdf(paste("Module_", input$selectModuleSec ,"_Relative_Abundance.pdf", sep=""), width = input$widthPDF_D2, height = input$heightPDF_D2)
        
        print(ggplot(module_Relative_abundanceSec(), aes_string(x = "variable", y = "value", fill = input$selectTaxo2Sec)) +
                geom_bar(stat = "identity") +
                labs(x = "Samples", y = "Relative Abundance", title = paste("Relative abundance at the ", input$selectTaxo2Sec, " Taxonomic level of the ", input$selectModuleSec, " module.", sep = "")) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)))
        
        dev.off()
        fs <- c(fs, paste("Module_", input$selectModuleSec ,"_Relative_Abundance.pdf", sep=""))
      }
      
      
      pdf("Module_Correlation_External_Trait.pdf", width = input$widthPDF_D2, height = input$heightPDF_D2)
      print(ggplot(data = selectedCondDFSec(), aes(x = rn, y = variable, fill = value)) +
              geom_tile(color = "white") +
              scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                   midpoint = 0, limit = c(-1,1), space = "Lab",
                                   name=paste("Module\n Conditions \nCorrelation"))+
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                               size = 10, hjust = 1))+
              geom_text(aes(label = label), color = "black", size = 3) +
              coord_fixed() +
              xlab('Module') +
              ylab(''))
      dev.off()
      fs <- c(fs, "Module_Correlation_External_Trait.pdf")
    }else{
      for (i in 1:ncol(selectedMEs2())){
        svg(paste("Contribution_MEs_", substr(colnames(selectedMEs2())[i], 3, 40), ".svg", sep = ""), width = input$widthPDF_D2, height = input$heightPDF_D2)
        print(barplot(selectedMEs2()[,i], col= substr(colnames(selectedMEs2())[i], 3, 40), main="", cex.main=2,
                      ylab="eigengene expression", cex.names = 0.7, names.arg = rownames(sampleAnnot_sec()), las = 2))
        fs <- c(fs, paste("Contribution_MEs_", substr(colnames(selectedMEs2())[i], 3, 40), ".svg", sep = ""))
      }
      svg(paste("Module_", input$selectModuleSec ,"_Corr_Vs_Membership.svg", sep=""), width = input$widthPDF_D2, height = input$heightPDF_D2)
      print(ggplot(moduleMembershipSec(), aes(x = moduleMembership, y = TraitSignificance)) +
              geom_point(color = input$selectModuleSec) +
              geom_smooth(method=lm, color = "white") +
              labs(x = paste("Module Membership in", input$selectModuleSec, "module"), y = paste("OTU significance for", input$selectVariable1Sec),
                   title = paste("Module membership vs. gene significance\n", textModuleMembershipSec() ,sep = "")) +
              geom_density_2d(color = '#E69F00'))
      
      dev.off()
      fs <- c(fs, paste("Module_", input$selectModuleSec ,"_Corr_Vs_Membership.svg", sep=""))
      if (input$TaxonFile1){
        svg(paste("Module_", input$selectModuleSec ,"_Relative_Abundance.svg", sep=""), width = input$widthPDF_D2, height = input$heightPDF_D2)
        
        print(ggplot(module_Relative_abundanceSec(), aes_string(x = "variable", y = "value", fill = input$selectTaxo2Sec)) +
                geom_bar(stat = "identity") +
                labs(x = "Samples", y = "Relative Abundance", title = paste("Relative abundance at the ", input$selectTaxo2Sec, " Taxonomic level of the ", input$selectModuleSec, " module.", sep = "")) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)))
        
        dev.off()
        fs <- c(fs, paste("Module_", input$selectModuleSec ,"_Relative_Abundance.svg", sep=""))
      }
      
      svg("Module_Correlation_External_Trait.svg", width = input$widthPDF_D2, height = input$heightPDF_D2)
      print(ggplot(data = selectedCondDFSec(), aes(x = rn, y = variable, fill = value)) +
              geom_tile(color = "white") +
              scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                   midpoint = 0, limit = c(-1,1), space = "Lab",
                                   name=paste("Module\n Conditions \nCorrelation"))+
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                               size = 10, hjust = 1))+
              geom_text(aes(label = label), color = "black", size = 3) +
              coord_fixed() +
              xlab('Module') +
              ylab(''))
      dev.off()
      fs <- c(fs, "Module_Correlation_External_Trait.svg")
    }
    
    zip(zipfile=filename, files=fs)
    
  },
  contentType = "application/zip")




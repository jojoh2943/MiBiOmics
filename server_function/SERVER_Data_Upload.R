#########################################################
#### SERVER FUNCTIONS : DATA UPLOAD AND MANIPULATION ####
#########################################################
# This script gathers all the function needed for the first tabulation of MiBiOmics and the data
# manipulation across the application


#### P1: REACTIVE OBJECTS ####

# Data Upload

sampleAnnot <- reactive({
  
  validate(
    need(input$file2$datapath != "", "Please select a data set")
  )
  if (input$rownames2 == FALSE){
    essai <- try(read.csv(input$file2$datapath,
                          header = input$header2,
                          sep = input$sep2,
                          dec = input$dec2))
    validate(
      need(substr(essai,1,5) != "Error", "Your annotation table cannot be uploaded. Please try to select a different separator or decimal. If you are still getting this error, verify the format of your file (you need a .csv file).")
    )
    sampleAnnot <- read.csv(input$file2$datapath,
                            header = input$header2,
                            sep = input$sep2,
                            dec = input$dec2)
    
  }else{
    essai <- try(read.csv(input$file2$datapath,
                          header = input$header2,
                          row.names = 1,
                          sep = input$sep2,
                          dec = input$dec2))
    validate(
      need(substr(essai,1,5) != "Error", "Your annotation table cannot be uploaded. Please try to select a different separator or decimal. If you are still getting this error, verify the format of your file (you need a .csv file).")
    )
    sampleAnnot <- read.csv(input$file2$datapath,
                            header = input$header2,
                            row.names = 1,
                            sep = input$sep2,
                            dec = input$dec2)
  }
  sampleAnnot
})

sampleAnnot_Default <- reactive({
  validate(
    need((input$LoadExample == "Yes" | input$LoadExample2 =="Yes"), "")
  )
  sampleMetadata<-"data/DS_Food_intake.csv" #sample sheet path
  sampleAnnot <- read.csv(sampleMetadata, sep = ";", header = TRUE, row.names = 1, dec = ",")
  # sampleAnnot <- sampleAnnot[which(rownames(sampleAnnot) %in% rownames(exprDat_Default())),]
  # sampleAnnot <- sampleAnnot[which(rownames(sampleAnnot) %in% rownames(exprDatSec_Default())),]
  #sampleAnnot <- sampleAnnot[which(rownames(sampleAnnot) %!in% input$selectSample),]
  sampleAnnot
})

sampleName <- reactive({
  if (input$LoadExample == "Yes" | input$LoadExample2 == "Yes" ){
    sampleName <- rownames(sampleAnnot_Default())
  }else{
    sampleName <- rownames(sampleAnnot())
  }
  sampleName
})

sampleAnnot1 <- reactive({
  if (input$LoadExample == "Yes" | input$LoadExample2 == "Yes" ){
    sampleAnnot <- sampleAnnot_Default()[which(rownames(sampleAnnot_Default()) %!in% input$selectSample),]
  }else{
    sampleAnnot <- sampleAnnot()[which(rownames(sampleAnnot()) %!in% input$selectSample),]
  }
})

sampleAnnot_Batch <- reactive({
  if (input$LoadExample == "Yes" | input$LoadExample2 == "Yes" ){
    sampleAnnot <- sampleAnnot_Default()[which(rownames(sampleAnnot_Default()) %!in% input$selectSample),]
  }else{
    sampleAnnot <- sampleAnnot()[which(rownames(sampleAnnot()) %!in% input$selectSample),]
  }
})

# Data Upload

exprDat <- reactive({
  validate(
    need(input$file1$datapath != "" && (input$LoadExample == "No" | input$LoadExample2 == "No"), "Please select a data set")
  )
  
  ### IMPORT CSV FORMAT
  if (input$rownames == FALSE){
    essai <- try(read.csv(input$file1$datapath,
                          header = input$header,
                          sep = input$sep,
                          dec = input$dec))
    validate(
      need(substr(essai,1,5) != "Error", "Your counting table cannot be uploaded. Please try to select a different separator or decimal. If you are still getting this error, verify the format of your file (you need a .csv file).")
    )
    exprDat <- read.csv(input$file1$datapath,
                        header = input$header,
                        sep = input$sep,
                        dec = input$dec)
  }else{
    essai <- try(read.csv(input$file1$datapath,
                          header = input$header,
                          row.names = 1,
                          sep = input$sep,
                          dec = input$dec))
    validate(
      need(substr(essai,1,5) != "Error", "Your counting table cannot be uploaded. Please try to select a different separator or decimal. If you are still getting this error, verify the format of your file (you need a .csv file).")
    )
    exprDat <- read.csv(input$file1$datapath,
                        header = input$header,
                        row.names = 1,
                        sep = input$sep,
                        dec = input$dec)
  }
  # Look if rows are OTUs and columns are samples.
  if (length(exprDat[which(colnames(exprDat) %in% rownames(sampleAnnot()))]) > 2){
    exprDat <- as.data.frame(t(exprDat))
  }
  exprDat[is.na(exprDat)] <- 0
  sampleID <- rownames(exprDat)
  # If the values are not numeric --> tranform values into numbers.
  if (!is.numeric(exprDat[,1])){
    exprDat <- sapply(exprDat, as.numeric)
    showNotification("The counting table does not contains numeric values. The values were transformed into numbers but your analysis may not worked as expected. Please check the format of your counting table.", type = "warning", duration = NULL)
  }
  if(ncol(exprDat) > 1000){
    showNotification(paste("The counting table contains", ncol(exprDat)," genes/OTUs. This might affect the speed of the analysis. You won't be able to see the file interactively but you will be able to download all the results."), type = "warning", duration = NULL)    
  }
  rownames(exprDat) <- sampleID
  
  
  
  
  ### SELECTED FILTRATION
  if (input$Filtration == "Yes"){
    if (input$TypeFiltration == "prevalence"){
      if(input$prevalence != 0){
        Prevalence_threshold = round((input$prevalence/100) * nrow(exprDat))
        names_OTU <- c()
        for (col in 1:ncol(exprDat)){
          nb_zero = 0
          for (row in 1:nrow(exprDat)){
            if (exprDat[row, col] == 0){
              nb_zero = nb_zero + 1
            }
          }
          
          if (nb_zero > (nrow(exprDat) - Prevalence_threshold)){
            names_OTU <- c(names_OTU, colnames(exprDat)[col])
          }
        }
        exprDat <- exprDat[,which(colnames(exprDat) %!in% names_OTU)]
      }
    }else{
      min_count = input$count
      exprDat <- exprDat[, which(colSums(exprDat) >= min_count)]
    }
  }
  
  # DATA NORMALIZATION
  if(input$Normalisation == "Yes"){
    if (input$TypeNormalisation == 'CSS'){
      MR_exprDat <- newMRexperiment(exprDat)
      p = cumNormStat(MR_exprDat) #default is 0.5
      MR_exprDat = cumNorm(MR_exprDat, p=p)
      exprDat <- MRcounts(MR_exprDat, norm = TRUE)
    }else{
      if (input$TypeNormalisation == 'TSS'){
        #data.TSS = t(apply(data.filter, 1, TSS.divide))
        exprDat <- apply(exprDat, 2, TSS.divide)
      }
    }
  }
  exprDat <- exprDat[which(rownames(exprDat) %in% rownames(sampleAnnot())),]
  exprDat <- exprDat[which(rownames(exprDat) %!in% input$selectSample),]
  sampleAnnotation <- sampleAnnot()
  sampleAnnotation <- sampleAnnotation[which(rownames(sampleAnnotation) %!in% input$selectSample),]
  # DATA TRANFORMATION
  if (input$Transformation == "Yes"){
    if (input$TypeTransformation == "Log10"){
      exprDat <- log((1 + exprDat), base = 10)
    }else{
      if(input$TypeTransformation == "Log2"){
        exprDat <- log((1 + exprDat), base = 2)
      }else{
        if(input$TypeTransformation == "Square"){
          exprDat <- exprDat^2
        }else{
          if(input$TypeTransformation == "sqrt"){
            exprDat <- sqrt(exprDat)
          }else{
            if(input$TypeTransformation == "Hellinger"){
              exprDat <- decostand(exprDat, method = "hellinger")
            }else{
              if (input$TypeTransformation == "CLR"){
                exprDat = clr(exprDat)
              }else{
                if (input$TypeTransformation == "ILR"){
                  exprDat.ilr = ilr(exprDat)
                  B <- exp(ilrBase(D=ncol(exprDat)))
                  exprDat = t(apply(exprDat.ilr, 1, function(x){
                    B[,1]^x[1] * B[,2]^x[2]
                  }))
                }
              }
            }
          }
        }
      }
    }
  }
  exprDat <- exprDat[which(rownames(exprDat) %in% rownames(sampleAnnotation)),]
  sampleAnnotation <- sampleAnnotation[which(rownames(sampleAnnotation) %in% rownames(exprDat)),]
  exprDat <- exprDat[match(rownames(sampleAnnotation), rownames(exprDat)),]
  # Batch effect 
  if (input$batchEffect){
    batch <- sampleAnnotation[,input$batchCol]
    
    # Create a design matrix
    modCombat <- model.matrix(~1, data=sampleAnnotation)
    #Apply ComBat function using parametric empirical Bayesian adjustments : returns a new expression matrix adjusted 
    #for batch
    exprDat <- ComBat(dat = t(as.matrix(exprDat)), batch = batch, mod=modCombat, par.prior = TRUE, prior.plots = FALSE)
    
    exprDat <- t(exprDat)
  }
  exprDat
})

exprDat_Default <- reactive({
  validate(
    need((input$LoadExample == 'Yes' | input$LoadExample2 =="Yes"), "")
  )
  expressionData<-"data/otu_table_silva.csv" #expression file path
  exprDat <- read.csv(expressionData, sep = ",", header = TRUE, row.names = 1)
  exprDat[is.na(exprDat)] <- 0
  ### SELECTED FILTRATION
  if (input$Filtration == "Yes"){
    if (input$TypeFiltration == "prevalence"){
      if(input$prevalence != 0){
        Prevalence_threshold = round((input$prevalence/100) * nrow(exprDat))
        names_OTU <- c()
        for (col in 1:ncol(exprDat)){
          nb_zero = 0
          for (row in 1:nrow(exprDat)){
            if (exprDat[row, col] == 0){
              nb_zero = nb_zero + 1
            }
          }
          
          if (nb_zero > (nrow(exprDat) - Prevalence_threshold)){
            names_OTU <- c(names_OTU, colnames(exprDat)[col])
          }
        }
        exprDat <- exprDat[,which(colnames(exprDat) %!in% names_OTU)]
      }
    }else{
      min_count = input$count
      exprDat <- exprDat[, which(colSums(exprDat) >= min_count)]
    }
  }
  exprDat <- exprDat[which(rownames(exprDat) %!in% input$selectSample),]
  sampleAnnotation <- sampleAnnot_Default()
  sampleAnnotation <- sampleAnnotation[which(rownames(sampleAnnotation) %!in% input$selectSample),]
  
  # DATA NORMALIZATION
  if(input$Normalisation == "Yes"){
    if (input$TypeNormalisation == 'CSS'){
      MR_exprDat <- newMRexperiment(exprDat)
      p = cumNormStat(MR_exprDat) #default is 0.5
      MR_exprDat = cumNorm(MR_exprDat, p=p)
      exprDat <- MRcounts(MR_exprDat, norm = TRUE)
    }else{
      if (input$TypeNormalisation == 'TSS'){
        exprDat <- apply(exprDat, 2, TSS.divide)
      }
    }
  }
  
  # DATA TRANSFORMATION
  if (input$Transformation == "Yes"){
    if (input$TypeTransformation == "Log10"){
      exprDat <- log((1 + exprDat), base = 10)
    }else{
      if(input$TypeTransformation == "Log2"){
        exprDat <- log((1 + exprDat), base = 2)
      }else{
        if(input$TypeTransformation == "Square"){
          exprDat <- exprDat^2
        }else{
          if(input$TypeTransformation == "sqrt"){
            exprDat <- sqrt(exprDat)
          }else{
            if(input$TypeTransformation == "Hellinger"){
              exprDat <- decostand(exprDat, method = "hellinger")
            }else{
              if (input$TypeTransformation == "CLR"){
                exprDat = clr(exprDat)
                #print(exprDat)
              }else{
                if (input$TypeTransformation == "ILR"){
                  exprDat.ilr = ilr(exprDat)
                  B <- exp(ilrBase(D=ncol(exprDat)))
                  exprDat = t(apply(exprDat.ilr, 1, function(x){
                    B[,1]^x[1] * B[,2]^x[2]
                  }))
                }
              }
            }
          }
        }
      }
    }
  }
  exprDat <- exprDat[which(rownames(exprDat) %in% rownames(sampleAnnotation)),]
  sampleAnnotation <- sampleAnnotation[which(rownames(sampleAnnotation) %in% rownames(exprDat)),]
  exprDat <- exprDat[match(rownames(sampleAnnotation), rownames(exprDat)),]
  # Batch effect 
  if (input$batchEffect){
    batch <- sampleAnnotation[,input$batchCol]
    
    # Create a design matrix
    modCombat <- model.matrix(~1, data=sampleAnnotation)
    #Apply ComBat function using parametric empirical Bayesian adjustments : returns a new expression matrix adjusted 
    #for batch
    
    exprDat <- ComBat(dat = t(as.matrix(exprDat)), batch = batch, mod=modCombat,par.prior = FALSE, prior.plots = FALSE)
    exprDat <- t(exprDat)
  }
  
  exprDat
})

exprDat_present <- reactive({
  if (input$LoadExample == "Yes" | input$LoadExample2 == "Yes"){
    exprDat_present <- exprDat_Default()
  }else{
    validate(
      need(input$file1$datapath != "", "Please select a data set")
    )
    exprDat_present <- exprDat()     
  }
  if (nchar(colnames(exprDat_present)[5]) > 20 ){
    for (i in 1:ncol(exprDat_present)){
      colnames(exprDat_present)[i] <- paste(input$CountingT, i, sep = "") 
    }
  }
  exprDat_present
})


exprDatSec <- reactive ({
  validate(
    need(input$file4$datapath != "", "Please select a data set")
  )
  ### IMPORT CSV FORMAT
  if (input$rownames4 == FALSE){
    essai <- try(read.csv(input$file4$datapath,
                          header = input$header4,
                          sep = input$sep4,
                          dec = input$dec4))
    validate(
      need(substr(essai,1,5) != "Error", "Your second omic table cannot be uploaded. Please try to select a different separator or decimal. If you are still getting this error, verify the format of your file (you need a .csv file).")
    )
    exprDatSec <- read.csv(input$file4$datapath,
                           header = input$header4,
                           sep = input$sep4,
                           dec = input$dec4)
  }else{
    essai <- try(read.csv(input$file4$datapath,
                          header = input$header4,
                          row.names = 1,
                          sep = input$sep4,
                          dec = input$dec4))
    validate(
      need(substr(essai,1,5) != "Error", "Your second omic table cannot be uploaded. Please try to select a different separator or decimal. If you are still getting this error, verify the format of your file (you need a .csv file).")
    )
    exprDatSec <- read.csv(input$file4$datapath,
                           header = input$header4,
                           row.names = 1,
                           sep = input$sep4,
                           dec = input$dec4)
  }
  if (length(exprDatSec[which(colnames(exprDatSec) %in% rownames(sampleAnnot()))]) > 2){
    exprDatSec <- as.data.frame(t(exprDatSec))
  }
  exprDatSec <- exprDatSec[which(rownames(exprDatSec) %in% rownames(sampleAnnot())),]
  validate(
    need(nrow(exprDatSec) != 0, "The sample names are differents between the datasets.")
  )
  
  exprDatSec[is.na(exprDatSec)] <- 0
  sampleID <- rownames(exprDatSec)
  # If the values are not numeric --> tranform values into numbers.
  if (!is.numeric(exprDatSec[,1])){
    exprDatSec <- sapply(exprDatSec, as.numeric)
    showNotification("The counting table does not contains numeric values. The values were transformed into numbers but your analysis may not worked as expected. Please check the format of your counting table.", type = "warning", duration = NULL)
  }else{
    exprDatSec <- sapply(exprDatSec, as.numeric)
  }
  if(ncol(exprDatSec) > 1000){
    showNotification(paste("The counting table contains", ncol(exprDatSec)," genes/OTUs. This might affect the speed of the analysis. You won't be able to see the file interactively but you will be able to download all the results."), type = "warning", duration = NULL)    
  }
  rownames(exprDatSec) <- sampleID
  
  
  ### SELECTED FILTRATION
  if (input$Filtration1 == "Yes"){
    if (input$TypeFiltration1 == "prevalence"){
      if(input$prevalence1 != 0){
        Prevalence_threshold = round((input$prevalence1/100) * nrow(exprDatSec))
        names_OTU <- c()
        for (col in 1:ncol(exprDatSec)){
          nb_zero = 0
          for (row in 1:nrow(exprDatSec)){
            if (exprDatSec[row, col] == 0){
              nb_zero = nb_zero + 1
            }
          }
          
          if (nb_zero > (nrow(exprDatSec) - Prevalence_threshold)){
            names_OTU <- c(names_OTU, colnames(exprDatSec)[col])
          }
        }
        exprDatSec <- exprDatSec[,which(colnames(exprDatSec) %!in% names_OTU)]
      }
    }else{
      min_count = input$count1
      exprDatSec <- exprDatSec[, which(colSums(exprDatSec) >= min_count)]
    }
  }
  
  # DATA NORMALIZATION
  if(input$Normalisation1 == "Yes"){
    if (input$TypeNormalisation1 == 'CSS'){
      MR_exprDatSec <- newMRexperiment(exprDatSec)
      if (input$LoadExample2 == "Yes" || input$LoadExample == "Yes" || input$OmicTable != 'OTUs'){
        p <- 0.5
      }else{
        p = cumNormStat(MR_exprDatSec) #default is 0.5
      }
      MR_exprDatSec = cumNorm(MR_exprDatSec, p=p)
      exprDatSec <- MRcounts(MR_exprDatSec, norm = TRUE)
    }else{
      if (input$TypeNormalisation1 == 'TSS'){
        exprDatSec <- apply(exprDatSec, 2, TSS.divide)
      }
    }
  }
  exprDatSec <- exprDatSec[which(rownames(exprDatSec) %in% rownames(sampleAnnot())),]
  exprDatSec <- exprDatSec[which(rownames(exprDatSec) %!in% input$selectSample),]
  sampleAnnotation <- sampleAnnot()
  sampleAnnotation <- sampleAnnotation[which(rownames(sampleAnnotation) %!in% input$selectSample),]
  
  # DATA TRANFORMATION
  if (input$Transformation1 == "Yes"){
    if (input$TypeTransformation1 == "Log10"){
      exprDatSec <- log((1 + exprDatSec), base = 10)
    }else{
      if(input$TypeTransformation1 == "Log2"){
        exprDatSec <- log((1 + exprDatSec), base = 2)
      }else{
        if(input$TypeTransformation1 == "Square"){
          exprDatSec <- exprDatSec^2
        }else{
          if(input$TypeTransformation1 == "sqrt"){
            exprDatSec <- sqrt(exprDatSec)
          }else{
            if(input$TypeTransformation1 == "Hellinger"){
              exprDatSec <- decostand(exprDatSec, method = "hellinger")
            }else{
              if (input$TypeTransformation1 == "CLR"){
                exprDatSec = clr(exprDatSec)
              }else{
                if (input$TypeTransformation1 == "ILR"){
                  exprDatSec.ilr = ilr(exprDatSec)
                  B <- exp(ilrBase(D=ncol(exprDatSec)))
                  exprDatSec = t(apply(exprDatSec.ilr, 1, function(x){
                    B[,1]^x[1] * B[,2]^x[2]
                  }))
                }
              }
            }
          }
        }
      }
    }
  }
  exprDatSec[is.na(exprDatSec)] <- 0
  sampleAnnotation <- sampleAnnotation[which(rownames(sampleAnnotation) %in% rownames(exprDatSec)),]
  exprDatSec <- exprDatSec[which(rownames(exprDatSec) %in% rownames(sampleAnnotation)),]
  exprDatSec <- exprDatSec[match(rownames(sampleAnnotation), rownames(exprDatSec)),]
  
  # Batch effect 
  if (input$batchEffect){
    batch <- sampleAnnotation[,input$batchCol]
    
    # Create a design matrix
    modCombat <- model.matrix(~1, data=sampleAnnotation)
    #Apply ComBat function using parametric empirical Bayesian adjustments : returns a new expression matrix adjusted 
    #for batch
    
    exprDatSec <- ComBat(dat = t(as.matrix(exprDatSec)), batch = batch, mod=modCombat,par.prior = FALSE, prior.plots = FALSE)
    exprDatSec <- t(exprDatSec)
  }
  exprDatSec
})

exprDatSec_Default <- reactive({
  validate(
    need(input$LoadExample2 =="Yes", "")
  )
  
  expressionData<-"data/metabolites.txt" #expression file path
  exprDatSec <- read.csv(expressionData, sep = ";", header = TRUE, row.names = 1)
  #exprDatSec <- exprDatSec[which(rownames(exprDatSec) %in% rownames(exprDat_present())),]
  validate(
    need(nrow(exprDatSec) != 0, "The sample names are differents between the datasets.")
  )
  exprDatSec[is.na(exprDatSec)] <- 0
  
  ### SELECTED FILTRATION
  if (input$Filtration1 == "Yes"){
    if (input$TypeFiltration1 == "prevalence"){
      if(input$prevalence1 != 0){
        Prevalence_threshold = round((input$prevalence1/100) * nrow(exprDatSec))
        names_OTU <- c()
        for (col in 1:ncol(exprDatSec)){
          nb_zero = 0
          for (row in 1:nrow(exprDatSec)){
            if (exprDatSec[row, col] == 0){
              nb_zero = nb_zero + 1
            }
          }
          
          if (nb_zero > (nrow(exprDatSec) - Prevalence_threshold)){
            names_OTU <- c(names_OTU, colnames(exprDatSec)[col])
          }
        }
        exprDatSec <- exprDatSec[,which(colnames(exprDatSec) %!in% names_OTU)]
      }
    }else{
      min_count = input$count1
      exprDatSec <- exprDatSec[, which(colSums(exprDatSec) >= min_count)]
    }
  }
  exprDatSec <- exprDatSec[which(rownames(exprDatSec) %!in% input$selectSample),]
  sampleAnnotation <- sampleAnnot_Default()
  sampleAnnotation <- sampleAnnotation[which(rownames(sampleAnnotation) %!in% input$selectSample),]
  sampleAnnotation <- sampleAnnotation[which(rownames(sampleAnnotation) %in% rownames(exprDatSec)),]
  exprDatSec <- exprDatSec[which(rownames(exprDatSec) %in% rownames(sampleAnnotation)),]
  
  # DATA NORMALIZATION
  if(input$Normalisation1 == "Yes"){
    if (input$TypeNormalisation1 == 'CSS'){
      MR_exprDatSec <- newMRexperiment(exprDatSec)
      if (input$LoadExample2 == "Yes" || input$LoadExample == "Yes" || input$OmicTable != 'OTUs'){
        p <- 0.5
      }else{
        p = cumNormStat(MR_exprDatSec) #default is 0.5
      }
      MR_exprDatSec = cumNorm(MR_exprDatSec, p=p)
      exprDatSec <- MRcounts(MR_exprDatSec, norm = TRUE)
      
    }else{
      if (input$TypeNormalisation1 == 'TSS'){
        exprDatSec <- apply(exprDatSec, 2, TSS.divide)
      }
    }
  }
  
  # DATA TRANSFORMATION
  if (input$Transformation1 == "Yes"){
    if (input$TypeTransformation1 == "Log10"){
      exprDatSec <- log((1 + exprDatSec), base = 10)
    }else{
      if(input$TypeTransformation1 == "Log2"){
        exprDatSec <- log((1 + exprDatSec), base = 2)
      }else{
        if(input$TypeTransformation1 == "Square"){
          exprDatSec <- exprDatSec^2
        }else{
          if(input$TypeTransformation1 == "sqrt"){
            exprDatSec <- sqrt(exprDatSec)
          }else{
            if(input$TypeTransformation1 == "Hellinger"){
              exprDatSec <- decostand(exprDatSec, method = "hellinger")
            }else{
              if (input$TypeTransformation1 == "CLR"){
                exprDatSec = clr(exprDatSec)
              }else{
                if (input$TypeTransformation1 == "ILR"){
                  exprDatSec.ilr = ilr(exprDatSec)
                  B <- exp(ilrBase(D=ncol(exprDatSec)))
                  exprDatSec = t(apply(exprDatSec.ilr, 1, function(x){
                    B[,1]^x[1] * B[,2]^x[2]
                  }))
                }
              }
            }
          }
        }
      }
    }
  }
  exprDatSec <- exprDatSec[match(rownames(sampleAnnotation), rownames(exprDatSec)),]
  # Batch effect 
  if (input$batchEffect){
    batch <- sampleAnnotation[,input$batchCol]
    
    # Create a design matrix
    modCombat <- model.matrix(~1, data=sampleAnnotation)
    #Apply ComBat function using parametric empirical Bayesian adjustments : returns a new expression matrix adjusted 
    #for batch
    
    exprDatSec <- ComBat(dat = t(as.matrix(exprDatSec)), batch = batch, mod=modCombat,par.prior = FALSE, prior.plots = FALSE)
    exprDatSec <- t(exprDatSec)
  }
  exprDatSec
})

exprDatSec_present <- reactive({
  if (input$LoadExample2 == "Yes"){
    exprDatSec_present <- exprDatSec_Default()
  }else{
    validate(
      need(input$file1$datapath != "", "Please select a data set")
    )
    exprDatSec_present <- exprDatSec()     
  }
  if (nchar(colnames(exprDatSec_present)[5]) > 20 ){
    for (i in 1:ncol(exprDatSec_present)){
      colnames(exprDatSec_present)[i] <- paste(input$OmicTable, i, sep = "") 
    }
  }
  exprDatSec_present
})

exprDatSec_2 <- reactive({
  if (input$LoadExample2 == "Yes"){
    exprDatSec_2 <- exprDatSec_Default()
  }else{
    validate(
      need(input$file4$datapath != "", "Please select a data set")
    )
    exprDatSec_2 <- exprDatSec()     
  }
  exprDatSec_2
})


taxTable <- reactive({
  if (input$TypeAnalysis == "simple"){
    validate(
      need(input$file3$datapath != "", "Please select a data set")
    )
    
    if (input$TaxonFile){
      if (input$rownames3 == FALSE){
        essai <- try(read.csv(input$file3$datapath,
                              header = input$header3,
                              sep = input$sep3))
        validate(
          need(substr(essai,1,5) != "Error", "Your counting table cannot be uploaded. Please try to select a different separator or decimal. If you are still getting this error, verify the format of your file (you need a .csv file).")
        )
        taxTable <- read.csv(input$file3$datapath,
                             header = input$header3,
                             sep = input$sep3)
      }else{
        essai <- try(read.csv(input$file3$datapath,
                              header = input$header3,
                              sep = input$sep3,
                              row.names = 1))
        validate(
          need(substr(essai,1,5) != "Error", "Your counting table cannot be uploaded. Please try to select a different separator or decimal. If you are still getting this error, verify the format of your file (you need a .csv file).")
        )
        taxTable <- read.csv(input$file3$datapath,
                             header = input$header3,
                             sep = input$sep3,
                             row.names = 1)    
      }
    }else{
      taxTable <- 0
    }
  }else{
    validate(
      need(input$file4$datapath != "", "Please select a data set")
    )
    
    if (input$TaxonFile1){
      if (input$rownames3 == FALSE){
        essai <- try(read.csv(input$file3$datapath,
                              header = input$header3,
                              sep = input$sep3))
        validate(
          need(substr(essai,1,5) != "Error", "Your counting table cannot be uploaded. Please try to select a different separator or decimal. If you are still getting this error, verify the format of your file (you need a .csv file).")
        )
        taxTable <- read.csv(input$file3$datapath,
                             header = input$header3,
                             sep = input$sep3)
      }else{
        essai <- try(read.csv(input$file3$datapath,
                              header = input$header3,
                              sep = input$sep3,
                              row.names = 1))
        validate(
          need(substr(essai,1,5) != "Error", "Your counting table cannot be uploaded. Please try to select a different separator or decimal. If you are still getting this error, verify the format of your file (you need a .csv file).")
        )
        taxTable <- read.csv(input$file3$datapath,
                             header = input$header3,
                             sep = input$sep3,
                             row.names = 1)    
      }
    }else{
      taxTable <- 0
    }
  }
  exprDat <- exprDat()
  taxTable <- taxTable[which(rownames(taxTable) %in% colnames(exprDat)),]
  taxTable <- taxTable[match(rownames(taxTable), colnames(exprDat)),]
  taxTable
  
})

taxTable_Default <- reactive({
  validate(
    need((input$LoadExample == "Yes" | input$LoadExample2 =="Yes"), "")
  )
  OTU_taxa_correspondance <- "data/otu_taxa_correspondance_silva.csv" # A correspondance table containing the OTU and its correspondant taxon
  OTU_Taxa_table <- read.csv(OTU_taxa_correspondance, sep = ",", header = TRUE, row.names = 1)
  exprDat <- exprDat_Default()
  OTU_Taxa_table <- OTU_Taxa_table[which(rownames(OTU_Taxa_table) %in% colnames(exprDat)),]
  OTU_Taxa_table <- OTU_Taxa_table[match(rownames(OTU_Taxa_table), colnames(exprDat)),]
  OTU_Taxa_table
})

taxTable_present <- reactive({
  if (input$LoadExample == "Yes" | input$LoadExample2 =="Yes"){
    taxTable_present <- taxTable_Default()
  }else{
    validate(
      need(input$file3$datapath != "", "Please select a data set")
    )
    if (taxTable() != 0){
      taxTable_present <- taxTable()
    }else{
      taxTable_present <- 0
    }
  }
  if (nchar(rownames(taxTable_present)) > 20){
    for (i in 1:nrow(taxTable_present)){
      rownames(taxTable_present)[i] <- paste("OTU", i, sep = "") 
    }      
  }
  
  taxTable_present
})

taxTable_rownames <- reactive({
  if (input$LoadExample == "Yes" || input$LoadExample2 == "Yes"){
    taxTable_row <- rownames(taxTable_Default())
  }else{
    taxTable_row <- rownames(taxTable())
  }
  taxTable_row
})

taxTable_report <- reactive({
  if (input$LoadExample == "Yes" || input$LoadExample2 == "Yes"){
    taxTable <- taxTable_Default()
  }else{
    taxTable <- taxTable()
  }
  
  taxTable
})



exprDat_report <- reactive({
  if (input$LoadExample == "Yes" || input$LoadExample2 == "Yes"){
    exprDat <- exprDat_Default()
  }else{
    exprDat <- exprDat()
  }
  
  exprDat
})

#### P2: REACTIVE OBJECTS ####



exprDat_2 <- reactive({
  if (input$LoadExample == "No" && input$LoadExample2 =="No"){
    validate(
      need(input$file1$datapath != "", "Please upload a counting table (it can be genes, OTUs or metabolites) in the general parameter tab"),
      need(input$file2$datapath != "", "Please upload a annotation file in the general parameter tab"),
      if (input$TaxonFile){
        need(input$file3$datapath != "", "Please upload the corresponding taxa file for the OTUs counting table file in the general parameter tab")
      }
    )
    expr <- exprDat()
    annot <- sampleAnnot()
  }else{
    expr <- exprDat_Default()
    annot <- sampleAnnot_Default()
  }

    exprDat1 <- expr
  
  if (any(rownames(annot) %!in% rownames(expr)) == TRUE ){
    sampleAnnot1 <- annot[which(rownames(annot ) %in% rownames(expr ) ),]
  }else{
    sampleAnnot1 <- annot
  }
  if (any(rownames(exprDat1) %!in% rownames(sampleAnnot1)) == TRUE ){
    exprDat1[which(rownames(exprDat1 ) %in% rownames(sampleAnnot1 ) ),]
  }else{
    exprDat1
  }
  
  
})


sampleAnnot_2 <- reactive({
  if (input$LoadExample == "No" && input$LoadExample2 =="No"){
    validate(
      need(input$file1$datapath != "", "Please upload a counting table (it can be genes, OTUs or metabolites) in the general parameter tab"),
      need(input$file2$datapath != "", "Please upload a annotation file in the general parameter tab"),
      if (input$TaxonFile == "OTUs"){
        need(input$file3$datapath != "", "Please upload the corresponding taxa file for the OTUs counting table file in the general parameter tab")
      }
    )
    annot <- sampleAnnot()
  }else{
    annot <- sampleAnnot_Default()
  }
  if (any(rownames(exprDat_2()) %!in% rownames(annot)) == TRUE ){
    exprDat1 <- exprDat2()[which(rownames(exprDat_2() ) %in% rownames(annot ) ),]
  }else{
    exprDat1 <- exprDat_2()
  }
  if (any(rownames(annot) %!in% rownames(exprDat1)) == TRUE ){
    annot[which(rownames(annot ) %in% rownames(exprDat1 ) ),]
  }else{
    annot
  }
  
})

sampleAnnot_sec <- reactive({
  sampleAnnot <- sampleAnnot_2()[which(rownames(sampleAnnot_2()) %in% rownames(exprDatSec_2())),]
  sampleAnnot
})

exprDatSec_3 <- reactive({
  expr <- exprDatSec_2()
  expr <- expr[which(rownames(expr) %in% rownames(sampleAnnot_sec())),]
  expr <- expr[match(rownames(expr), rownames(sampleAnnot_sec())),]
  expr
})





taxTable1 <- reactive({
  if (input$LoadExample == 'No' && input$LoadExample2 =='No'){
    taxTable <- taxTable()
    taxTable <- taxTable[which(rownames(taxTable) %in% colnames(exprDat_2())),]
    taxTable <- taxTable[match(rownames(taxTable), colnames(exprDat_2())),]
    taxTable$rn <- colnames(exprDat_2())
    
  }else{
    taxTable <- taxTable_Default()
    taxTable <- taxTable[which(rownames(taxTable) %in% colnames(exprDat_2())),]
    taxTable <- taxTable[match(rownames(taxTable), colnames(exprDat_2())),]
    taxTable$rn <- colnames(exprDat_2())
    
  }
  taxTable
})

taxTable_Sec <- reactive({
  taxTable <- taxTable()
  taxTable <- taxTable[which(rownames(taxTable) %in% colnames(exprDatSec_3())),]
  taxTable <- taxTable[match(rownames(taxTable), colnames(exprDatSec_3())),]
  taxTable$rn <- colnames(exprDatSec_3())
  taxTable
})

#### P3: REACTIVE OBJECTS ####
exprDat_WGCNA <- reactive({
  expr <- exprDat_2()
  if (nrow(expr) > 14){
    gsg = goodSamplesGenes(expr, verbose = 3)
    #If gsg$allOK is TRUE, all genes have passed the cuts.  If not, we remove the offending genes and samples

    if (!gsg$allOK)
    {
      # Optionally, print the gene and sample names that were removed:
      if (sum(!gsg$goodGenes)>0)
      {
        printFlush(paste("Removing genes:", paste(colnames(expr)[!gsg$goodGenes], collapse = ", ")));
      }
      if (sum(!gsg$goodSamples)>0)
      {
        printFlush(paste("Removing samples:", paste(rownames(expr)[!gsg$goodSamples], collapse = ", ")));
      }
      exprDat1 <- expr[rownames(expr)[gsg$goodSamples], colnames(expr)[gsg$goodGenes]]
    }else{
      exprDat1 <- expr
    }
  }else{
  exprDat1 <- expr
  }
  exprDat1
})

exprDatSec_WGCNA <- reactive({
  expr <- exprDatSec_3()
  if (nrow(expr) > 14){
    gsg = goodSamplesGenes(expr, verbose = 3)
    #If gsg$allOK is TRUE, all genes have passed the cuts.  If not, we remove the offending genes and samples

    if (!gsg$allOK)
    {
      # Optionally, print the gene and sample names that were removed:
      if (sum(!gsg$goodGenes)>0)
      {
        printFlush(paste("Removing genes:", paste(colnames(expr)[!gsg$goodGenes], collapse = ", ")));
      }
      if (sum(!gsg$goodSamples)>0)
      {
        printFlush(paste("Removing samples:", paste(rownames(expr)[!gsg$goodSamples], collapse = ", ")));
      }
      expr <- expr[rownames(expr)[gsg$goodSamples], colnames(expr)[gsg$goodGenes]]
    }
  }
  expr
})

#### P5: REACTIVE OBJECTS ####


exprDat_3 <- reactive({
  expr <- exprDat_2()
  expr2 <- exprDatSec_3()
  expr <- expr[which(rownames(expr) %in% rownames(expr2)),]
  expr2 <- expr2[which(rownames(expr2) %in% rownames(expr)),]
  expr <- expr[match(rownames(expr), rownames(expr2)),]
  if (input$CountingT1 == "OTUs" && input$OmicTable == "OTUs" && input$LoadExample2 == "No"){
    colnames(expr) <- paste("a", colnames(expr), sep ="")
  }
  expr
}) 

exprDatSec_4 <- reactive({
  expr <- exprDat_2()
  expr2 <- exprDatSec_3()
  expr <- expr[which(rownames(expr) %in% rownames(expr2)),]
  expr2 <- expr2[which(rownames(expr2) %in% rownames(expr)),]
  expr2 <- expr2[match(rownames(expr2), rownames(expr)),]
  if (input$CountingT1 == "OTUs" && input$OmicTable == "OTUs" && input$LoadExample2 == "No"){
    colnames(expr2) <- paste("a", colnames(expr2), sep ="")
  }
  expr2
}) 

taxTable2 <- reactive({
  taxTable <- taxTable1()
  if (input$CountingT1 == "OTUs" && input$OmicTable == "OTUs" && input$LoadExample2 == "No"){
    taxTable$rn <- paste("a", rownames(taxTable), sep ="")
  }
  taxTable
})

selectedMetab <- reactive({
  expr <- exprDatSec_4()
  
  
  if (input$selectModule2_p5 != "General"){
    modGenes = (selectedDynamicColor2()== input$selectModule2_p5)
    
    Module_metab <- expr[,modGenes]
    unique_col <- c()
    for (i in 1:ncol(Module_metab)){
      if (length(unique(Module_metab[,i])) == 1){ #if there is only one value for all the samples
        unique_col <- c(unique_col, colnames(Module_metab)[i])
      }
    }
    Module_metab <- Module_metab[, which(colnames(Module_metab) %!in% unique_col)]
  }else{
    Module_metab <- expr
    unique_col <- c()
    for (i in 1:ncol(Module_metab)){
      if (length(unique(Module_metab[,i])) == 1){ #if there is only one value for all the samples
        unique_col <- c(unique_col, colnames(Module_metab)[i])
      }
    }
    Module_metab <- Module_metab[, which(colnames(Module_metab) %!in% unique_col)]
  }
  Module_metab[which(is.na(Module_metab)),] <- 0
  
  Module_metab
})

sampleAnnot_3 <- reactive({
  annot <- sampleAnnot_2()[which(rownames(sampleAnnot_2()) %in% rownames(exprDatSec_4())),]
  annot <- annot[match(rownames(annot), rownames(exprDatSec_4())),]
  annot
})

selected_moduleOTU <- reactive({
  if (input$selectModule1_p5 != "General"){
    modGenes = (selectedDynamicColor()== input$selectModule1_p5)
    
    Module_OTU <- exprDat_3()[modGenes]
    unique_col <- c()
    for (i in 1:ncol(Module_OTU)){
      if (length(unique(Module_OTU[,i])) == 1){ #if there is only one value for all the samples
        unique_col <- c(unique_col, colnames(Module_OTU)[i])
      }
    }
    Module_OTU <- Module_OTU[, which(colnames(Module_OTU) %!in% unique_col)]
  }else{
    Module_OTU <- exprDat_3()
    unique_col <- c()
    for (i in 1:ncol(Module_OTU)){
      if (length(unique(Module_OTU[,i])) == 1){ #if there is only one value for all the samples
        unique_col <- c(unique_col, colnames(Module_OTU)[i])
      }
    }
    Module_OTU <- Module_OTU[, which(colnames(Module_OTU) %!in% unique_col)]
  } 
  Module_OTU
})

#### HEATMAP 

exprDat_4 <- reactive({
  if (input$selectModule11_p5 == "Individual_Variables"){
    validate(
      need(input$selectVariablesExprDat != "", "Select variables to plot in the heatmap")
    )
    exprDat_3()[,input$selectVariablesExprDat]
  }else{
    modGenes = (selectedDynamicColor()== input$selectModule11_p5)
    exprDat_3()[,modGenes]
  }
})

exprDatSec_5 <- reactive({
  if (input$selectModule22_p5 == "Individual_Variables"){
    validate(
      need(input$selectVariables != "", "Select variables to plot in the heatmap")
    )
    exprDatSec_4()[,input$selectVariables]
  }else{
    modGenes = (selectedDynamicColor2()== input$selectModule22_p5)
    exprDatSec_4()[,modGenes] 
  }
  
})

#### INTERFACE VARIABLES ####

output$ViewPanel <- renderUI({
  if (input$TypeAnalysis == 'simple' && (input$TaxonFile == TRUE | (input$LoadExample == "Yes" ))){
    radioButtons("View1",
                 "View: ",
                 choices = c("Counting Table" = "counting", "Annotation Table" = "annot", "Taxa Table" = "taxa"),
                 selected = character(0))
  }else{
    if (input$TypeAnalysis == 'multivariate' && (input$TaxonFile1 == TRUE | input$LoadExample2 == "Yes")){
      radioButtons("View3",
                   "View: ",
                   choices = c("Counting Table" = "counting", "Annotation Table" = "annot", "Taxa Table" = "taxa", "Second Omic Table" = "sec"),
                   selected = character(0))
    }else{
      if (input$TypeAnalysis == 'multivariate' && input$TaxonFile1 == FALSE){
        radioButtons("View4",
                     "View: ",
                     choices = c("Counting Table" = "counting", "Annotation Table" = "annot", "Second Omic Table" = "sec"),
                     selected = character(0))
      }else{
        radioButtons("View2", 
                     "View: ", 
                     choices = c("Counting Table" = "counting", "Annotation Table" = "annot"),
                     selected = character(0))
      }
    }
  }
})

#select sample to remove
output$selectSample <- renderUI({
  if (input$LoadExample == "Yes" | input$LoadExample2 == "Yes"){
    selectInput("selectSample", "Select samples to remove: ",
                choices = sampleName(),
                multiple = TRUE)
  }else{
    selectInput("selectSample", "Select samples to remove: ",
                choices = rownames(sampleAnnot()),
                multiple = TRUE) 
  }
  
})

output$batchCol <-  renderUI({
  selectInput("batchCol", 
              label = "Choose the batch Column: ", 
              choices = colnames(sampleAnnot_Batch()), 
              selected = colnames(sampleAnnot_Batch())[2])
})

#### FIGURES OUTPUTS ####

output$ViewTable <- renderDataTable({
  if(is.null(input$View1)&is.null(input$View2)&is.null(input$View3)&is.null(input$View4)){
    validate(
      need(input$file1$datapath != "", "Please upload a counting table")
    )
    validate(
      need(input$file2$datapath != "", "Please upload an annotation table")
    )
    if(input$TypeAnalysis == 'multivariate'){
      if(input$TaxonFile1){
        validate(
          need(input$file3$datapath != "", "Please upload a Taxon File")
        )
        validate(
          need(input$View3 != "", "Please select a viewing option")
        )
      }
      validate(
        need(input$file4$datapath != "", "Please upload a Second omic table file")
      )
      validate(
        need(input$View4 != "", "Please select a viewing option")
      )
      
    }else{
      if(input$TaxonFile){
        validate(
          need(input$file3$datapath != "", "Please upload a Taxon File")
        )
        validate(
          need(input$View2 != "", "Please select a viewing option")
        )
      }else{
        validate(
          need(input$View1 != "", "Please select a viewing option")
        ) 
      }   
    }
    
  }else{
    if(input$TypeAnalysis == 'simple'){
      if(input$TaxonFile | input$LoadExample == "Yes"){
        if (input$View1 == "counting"){
          if (ncol(exprDat_present()) < 1000){
            datatable(exprDat_present()[,1:10],
                      options = list(lengthMenu = FALSE, pageLength = 5,  dom = 'tip'),
                      class = 'cell-border stripe')
          }else{
            datatable(head(exprDat_present()[,1:10]),
                      options = list(lengthMenu = FALSE, pageLength = 6,  dom = 'tip'),
                      class = 'cell-border stripe')
          }
          
        }else{
          if (input$View1 == "annot"){
            datatable(sampleAnnot_2(), 
                      options = list(lengthMenu = FALSE, pageLength = 5,  dom = 'tip'),
                      class = 'cell-border stripe')
          }else{
            datatable(taxTable_present(),
                      options = list(lengthMenu = FALSE, pageLength = 5, dom = 'tip'), 
                      class = 'cell-border stripe')
          }
        }
      }else{
        if (input$View2 == "annot"){
          datatable(sampleAnnot_2(), 
                    options = list(lengthMenu = FALSE, pageLength = 5,  dom = 'tip'),
                    class = 'cell-border stripe')
        }else{
          if(ncol(exprDat_present()) < 1000){
            datatable(exprDat_present()[,1:10],
                      options = list(lengthMenu = FALSE, pageLength = 5,  dom = 'tip'),
                      class = 'cell-border stripe')       
          }else{
            datatable(head(exprDat_present()[,1:10]),
                      options = list(lengthMenu = FALSE, pageLength = 5,  dom = 'tip'),
                      class = 'cell-border stripe')                
          }
        }
      }
    }else{
      if(input$TaxonFile1 | input$LoadExample2 == "Yes"){
        if (input$View3 == "counting"){
          if (ncol(exprDat_present()) < 1000){
            datatable(exprDat_present()[,1:10],
                      options = list(lengthMenu = FALSE, pageLength = 5,  dom = 'tip'),
                      class = 'cell-border stripe')
          }else{
            datatable(head(exprDat_present()[,1:10]),
                      options = list(lengthMenu = FALSE, pageLength = 6,  dom = 'tip'),
                      class = 'cell-border stripe')
          }
          
        }else{
          if (input$View3 == "annot"){
            datatable(sampleAnnot_2(), 
                      options = list(lengthMenu = FALSE, pageLength = 5,  dom = 'tip'),
                      class = 'cell-border stripe')
          }else{
            if (input$View3 == "taxa"){
              datatable(taxTable_present(),
                        options = list(lengthMenu = FALSE, pageLength = 5, dom = 'tip'), 
                        class = 'cell-border stripe')
            }else{
              datatable(exprDatSec_present()[,1:10],
                        options = list(lengthMenu = FALSE, pageLength = 5, dom = 'tip'), 
                        class = 'cell-border stripe')
            }
          }
        }
      }else{
        if (input$View4 == "annot"){
          datatable(sampleAnnot_2(), 
                    options = list(lengthMenu = FALSE, pageLength = 5,  dom = 'tip'),
                    class = 'cell-border stripe')
        }else{
          if (input$View4 =="counting"){
            if(ncol(exprDat_present()) < 1000){
              datatable(exprDat_present()[,1:10],
                        options = list(lengthMenu = FALSE, pageLength = 5,  dom = 'tip'),
                        class = 'cell-border stripe')       
            }else{
              datatable(head(exprDat_present()[,1:10]),
                        options = list(lengthMenu = FALSE, pageLength = 5,  dom = 'tip'),
                        class = 'cell-border stripe')                
            }
          }else{
            datatable(exprDatSec_present()[,1:10],
                      options = list(lengthMenu = FALSE, pageLength = 5, dom = 'tip'), 
                      class = 'cell-border stripe')
          }
        }
      }
    }
  }
  
  
})

output$ViewDim <- renderDataTable({
  if(is.null(input$View1)&is.null(input$View2)&is.null(input$View3)&is.null(input$View4)){
    validate(
      need(input$file1$datapath != "", "Please upload a counting table")
    )
    validate(
      need(input$file2$datapath != "", "Please upload an annotation table")
    )
    if(input$TypeAnalysis == 'multivariate'){
      if(input$TaxonFile1){
        validate(
          need(input$file3$datapath != "", "Please upload a Taxon File")
        )
        validate(
          need(input$View3 != "", "Please select a viewing option")
        )
      }
      validate(
        need(input$file4$datapath != "", "Please upload a Second omic table file")
      )
      validate(
        need(input$View4 != "", "Please select a viewing option")
      )
      
    }else{
      if(input$TaxonFile){
        validate(
          need(input$file3$datapath != "", "Please upload a Taxon File")
        )
        validate(
          need(input$View2 != "", "Please select a viewing option")
        )
      }else{
        validate(
          need(input$View1 != "", "Please select a viewing option")
        ) 
      }   
    }
    
  }else{
    if(input$TypeAnalysis == 'simple'){
      if(input$TaxonFile | input$LoadExample == "Yes"){
        if (input$View1 == "counting"){
          dimension <- as.data.frame(dim(exprDat_present()))
          colnames(dimension) <- c("dimension")
          rownames(dimension) <- c("row", "col")
          datatable(dimension,
                    options = list(lengthMenu = FALSE, dom = 'tip'),
                    class = 'cell-border stripe')
          
        }else{
          if (input$View1 == "annot"){
            dimension <- as.data.frame(dim(sampleAnnot_2()))
            colnames(dimension) <- c("dimension")
            rownames(dimension) <- c("row", "col")
            datatable(dimension,
                      options = list(lengthMenu = FALSE, dom = 'tip'),
                      class = 'cell-border stripe')
          }else{
            dimension <- as.data.frame(dim(taxTable_present()))
            colnames(dimension) <- c("dimension")
            rownames(dimension) <- c("row", "col")
            datatable(dimension,
                      options = list(lengthMenu = FALSE, dom = 'tip'),
                      class = 'cell-border stripe')
          }
        }
      }else{
        if (input$View2 == "annot"){
          dimension <- as.data.frame(dim(sampleAnnot_2()))
          colnames(dimension) <- c("dimension")
          rownames(dimension) <- c("row", "col")
          datatable(dimension,
                    options = list(lengthMenu = FALSE, dom = 'tip'),
                    class = 'cell-border stripe')
        }else{
          dimension <- as.data.frame(dim(exprDat_present()))
          colnames(dimension) <- c("dimension")
          rownames(dimension) <- c("row", "col")
          datatable(dimension,
                    options = list(lengthMenu = FALSE, dom = 'tip'),
                    class = 'cell-border stripe')
        }
      }
    }else{
      if(input$TaxonFile1 | input$LoadExample2 == "Yes"){
        if (input$View3 == "counting"){
          dimension <- as.data.frame(dim(exprDat_present()))
          colnames(dimension) <- c("dimension")
          rownames(dimension) <- c("row", "col")
          datatable(dimension,
                    options = list(lengthMenu = FALSE, dom = 'tip'),
                    class = 'cell-border stripe')
          
        }else{
          if (input$View3 == "annot"){
            dimension <- as.data.frame(dim(sampleAnnot_2()))
            colnames(dimension) <- c("dimension")
            rownames(dimension) <- c("row", "col")
            datatable(dimension,
                      options = list(lengthMenu = FALSE, dom = 'tip'),
                      class = 'cell-border stripe')
          }else{
            if (input$View3 == "taxa"){
              dimension <- as.data.frame(dim(taxTable_present()))
              colnames(dimension) <- c("dimension")
              rownames(dimension) <- c("row", "col")
              datatable(dimension,
                        options = list(lengthMenu = FALSE, dom = 'tip'),
                        class = 'cell-border stripe')
            }else{
              dimension <- as.data.frame(dim(exprDatSec_present()))
              colnames(dimension) <- c("dimension")
              rownames(dimension) <- c("row", "col")
              datatable(dimension,
                        options = list(lengthMenu = FALSE, dom = 'tip'),
                        class = 'cell-border stripe')
            }
          }
        }
      }else{
        if (input$View4 == "annot"){
          dimension <- as.data.frame(dim(sampleAnnot_2()))
          colnames(dimension) <- c("dimension")
          rownames(dimension) <- c("row", "col")
          datatable(dimension,
                    options = list(lengthMenu = FALSE, dom = 'tip'),
                    class = 'cell-border stripe')
        }else{
          if (input$View4 =="counting"){
            dimension <- as.data.frame(dim(exprDat_present()))
            colnames(dimension) <- c("dimension")
            rownames(dimension) <- c("row", "col")
            datatable(dimension,
                      options = list(lengthMenu = FALSE, dom = 'tip'),
                      class = 'cell-border stripe')
          }else{
            dimension <- as.data.frame(dim(exprDatSec_present()))
            colnames(dimension) <- c("dimension")
            rownames(dimension) <- c("row", "col")
            datatable(dimension,
                      options = list(lengthMenu = FALSE, dom = 'tip'),
                      class = 'cell-border stripe')
          }
        }
      }
    }
  }
  
  
})
output$dim <- renderPrint({
  dim(exprDat_present())
})
output$dimAnnot <- renderPrint({
  dim(sampleAnnot_2())
})


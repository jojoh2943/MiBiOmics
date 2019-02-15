################################################
#### SERVER FUNCTIONS : HELP AND ABOUT PAGE ####
################################################


#### OBSERVE EVENTS ####





exampleDATA <- data.frame("OTU1" = c(0,0,1,6,0), "OTU2" = c(2,0,1,0,0), "OTU3" = c(7,0,4,6,0), "OTU4" = c(2,3,0,0,0), "OTU5" = c(7,0,0,0,0), "OTU6" = c(0,0,1,0,0))
rownames(exampleDATA) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
exampleANNOT <- data.frame("SampleSite" = c("Site1", "Site6", "Site1", "Site5", "Site3"), "AlphaDiversity" = c(0.56,0.72,0.11,0.32,0.06), "Depth" = c(40,32,56,106,12))
rownames(exampleANNOT) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
exampleTAXA <-  data.frame("Kingdom" = c("Bacteria","Bacteria","Bacteria","Bacteria","Bacteria", "Bacteria"), "Phylum" = c("Bacteroidetes","Bacteroidetes","Bacteroidetes","Proteobacteria","Proteobacteria", "Firmicutes"), "Class" = c("Bacteroidia","Bacteroidia","Bacteroidia","Deltaproteobacteria","Deltaproteobacteria", "Bacilli"), "Order" = c("Bacteroidales","Bacteroidales","Bacteroidales","Desulfovibrionales","Desulfovibrionales", "Lactobacillales"), "Family" = c("Porphyromonadaceae","Bacteroidaceae","Bacteroidaceae","Desulfovibrionaceae","Desulfovibrionaceae", "Lactobacillaceae"), "Genus" = c("Parabacteroides","Bacteroidaceae","Alistipes","Desulfovibrio","Desulfovibrio", "Lactobacillus"), "Species" =c("Unassigned", "Unassigned", "Uncultured bacterium", "Unassigned", "Unassigned", "Unassigned"))
rownames(exampleTAXA) <- c("OTU1", "OTU2", "OTU3", "OTU4", "OTU5", "OTU6")
#### HELP PAGE 
output$exampleTable <- DT::renderDataTable({
  DT::datatable(exampleDATA, options = list(lengthMenu = FALSE, dom = 'tip'),
                class = 'cell-border stripe',
                caption = 'Table 1: This is an example counting table.')
})

output$exampleAnnot <- DT::renderDataTable({
  DT::datatable(exampleANNOT, options = list(lengthMenu = FALSE, dom = 'tip'),
                class = 'cell-border stripe',
                caption = 'Table 2: This is an example annotation table.')
})

output$exampleTaxa <- DT::renderDataTable({
  DT::datatable(exampleTAXA, options = list(lengthMenu = FALSE, dom = 'tip'),
                class = 'cell-border stripe',
                caption = 'Table 3: This is an example taxonomic annotation table.')
})


output$sectionOne <- renderText({
  input$commentSection1
})

#### REPORT GENERATION ####

tempReport <- reactive({
  tempReport <- file.path(tempdir(), "report.Rmd")
  tempReport
})

generateReport <- observeEvent(input$produceReport,{
  input$produceReport
  withProgress({
    
    filename = "report.html"
    datasets <- list()
    if (input$LoadExample == "No" && input$LoadExample2 == "No"){
      datasets[["exprDat"]] <- exprDat_report()
      datasets[["sampleAnnot"]] <- sampleAnnot()
      if (input$TaxonFile){
        datasets[["taxaTable"]] <- taxTable()
        datasets[["type"]] <- "ST" # Simple + Taxon
      }else{
        datasets[["type"]] <- "S" # Simple without taxon
      }
      if (input$TypeAnalysis == "multivariate"){
        datasets[["exprDatSec"]] <- exprDatSec()
        if (input$TaxonFile1){
          datasets[["taxaTable"]] <- taxTable_report()
          datasets[["type"]] <- "MT" # Multivariate + Taxon
        }else{
          datasets[["type"]] <- "M" # Multivariate without taxon
        } 
      }
    }else{
      datasets[["exprDat"]] <- exprDat_report()
      datasets[["sampleAnnot"]] <- sampleAnnot_Default()
      datasets[["taxaTable"]] <- taxTable_report()
      datasets[["type"]] <- "ST" # Simple + Taxon 
      if (input$TypeAnalysis == "multivariate"){
        datasets[["exprDatSec"]] <- exprDatSec_Default()
        datasets[["type"]] <- "MT" # Simple + Taxon
      }
    }
    list_treatment <- list()
    #Create a vector containing all the inputs:
    #Input for filtration
    if (input$Filtration == "Yes"){
      list_treatment[["filtration"]] <- "F"
      if (input$TypeFiltration == "count"){
        list_treatment[["filtration"]] <- paste(list_treatment[["filtration"]], "C", sep = "")
        list_treatment[["filtration_value"]] <- input$count
        filtration_text <- paste("A filtration on the first counting table was realized based on count. A metabolite/OTU/gene is kept in the analysis if its total number of count thoughout the samples is superior to", input$count, ".")
      }else{
        list_treatment[["filtration"]] <- paste(list_treatment[["filtration"]], "P", sep = "")
        list_treatment[["filtration_value"]] <- input$prevalence
        filtration_text <- paste("A filtration on the first counting table was realized based on prevalence. A metabolite/OTU/gene is kept in the analysis if it is present in at least", input$prevalence, "% of the samples.")
      }
    }else{
      list_treatment[["filtration"]] <- "NONE"
      filtration_text <- ("")
    }
    if (datasets[["type"]] == "MT" || datasets[["type"]] == "M"){
      if (input$Filtration1 == "Yes"){
        list_treatment[["filtrationD2"]] <- "F"
        if (input$TypeFiltration1 == "count"){
          list_treatment[["filtrationD2"]] <- paste(list_treatment[["filtrationD2"]], "C", sep = "")
          list_treatment[["filtration_valueD2"]] <- input$count1
          filtration_text <- paste(filtration_text, "A filtration on the second counting table was realized based on count. A metabolite/OTU/gene is kept in the analysis if its total number of count thoughout the samples is superior to", input$count1)
        }else{
          list_treatment[["filtrationD2"]] <- paste(list_treatment[["filtrationD2"]], "P", sep = "")
          list_treatment[["filtration_valueD2"]] <- input$prevalence1
          filtration_text <- paste(filtration_text, "A filtration on the second counting table was realized based on prevalence. A metabolite/OTU/gene is kept in the analysis if it is present in at least", input$prevalence1, "% of the samples.")
        }  
      }else{
        list_treatment[["filtrationD2"]] <- "NONE"
        filtration_text <- paste(filtration_text,"")
      }
    }
    #Input for normalization
    if (input$Normalisation == "Yes"){
      if (input$TypeNormalisation == "MinMax"){
        list_treatment[["normalization"]] <- "MinMax"
        normalization_text <- paste("A normalization on the first counting table was realized based on the Min - Max method." )
      }else{
        list_treatment[["normalization"]] <- "Z"
        normalization_text <- paste("A normalization on the first counting table was realized based on the Z-score method (function scale() in R).")
      }
    }else{
      list_treatment[["normalization"]] <- "NONE"
      normalization_text <- ("")
    }
    if (datasets[["type"]] == "MT" || datasets[["type"]] == "M"){
      if (input$Normalisation1 == "Yes"){
        if (input$TypeNormalisation1 == "MinMax"){
          list_treatment[["normalizationD2"]] <- "MinMax"
          normalization_text <- paste(normalization_text, "A normalization on the second counting table was realized based on the Min - Max method.")
        }else{
          list_treatment[["normalizationD2"]] <- "Z"
          normalization_text <- paste(normalization_text, "A normalization on the second counting table was realized based on the Z-score method (function scale() in R).")
        }  
      }else{
        list_treatment[["normalizationD2"]] <- "NONE"
        normalization_text <- paste(normalization_text,"")
      }
    }
    
    #Input for transformation
    if (input$Transformation == "Yes"){
      list_treatment[["transformation"]] <- input$TypeTransformation
      transformation_text <- paste("A", input$TypeTransformation ,"tranformation on the first counting table was realized." )
    }else{
      list_treatment[["transformation"]] <- "NONE"
      transformation_text <- ("")
    }
    if (datasets[["type"]] == "MT" || datasets[["type"]] == "M"){
      if (input$Transformation1 == "Yes"){
        list_treatment[["transformationD2"]] <- input$TypeTransformation1
        transformation_text <- paste(transformation_text,"A", input$TypeTransformation1 ,"tranformation on the second counting table was realized." )
      }else{
        list_treatment[["transformationD2"]] <- "NONE"
        transformation_text <- paste(transformation_text,"")
      }
    }
    if (!is.null(input$selectSample)){
      list_treatment[["SampleToRemove"]] <- input$selectSample
      SampleToRemove <- paste(input$selectSample, collapse = ", ")
    }else{
      list_treatment[["SampleToRemove"]] <- "NONE"
      SampleToRemove <- "None"
    }
    vector_input <- c(filtration_text, 
                      normalization_text, 
                      transformation_text, 
                      SampleToRemove, 
                      input$commentSection1, 
                      input$LoadExample, 
                      input$LoadExample2, 
                      input$commentSection2,
                      input$Power_P6,
                      input$Power_D2_P6,
                      input$ModuleSize_P6,
                      input$ModuleSize_D2_P6,
                      input$signed_P6,
                      input$signed_D2_P6,
                      input$corNetwork_P6,
                      input$corNetwork_D2_P6,
                      input$commentSection3,
                      input$commentSection4,
                      input$commentSection5)
    incProgress(0.25, detail = "uploading parameters")
    filePath <- getwd()
    write_report(filePath, datasets, vector_input)
    incProgress(0.60, detail = "writing report")
    # Copy the report file to a temporary directory before processing it, in
    # case we don't have write permissions to the current working dir (which
    # can happen when deployed).
    if (input$TypeAnalysis == "simple"){
      list_input_P6 <- list(dist_PCoA_P2 = input$selectDist_P6, 
                            method_dendrogramme_P2 = input$selectMethod_P6,
                            power = input$Power_P6,
                            moduleSize = input$ModuleSize_P6,
                            signed = input$signed_P6,
                            cor_P3 = input$corNetwork_P6,
                            cor_P4 = input$selectCorrelation_P6,
                            loadExample = input$LoadExample,
                            typeData_D1 = input$CountingT)
    }else{
      list_input_P6 <- list(dist_PCoA_P2 = input$selectDist_P6, 
                            method_dendrogramme_P2 = input$selectMethod_P6, 
                            dist_PCoA_D2_P2 = input$selectDist_D2_P6,
                            method_dendrogramme_D2_P2 = input$selectMethod_D2_P6,
                            power = input$Power_P6,
                            power_D2 = input$Power_D2_P6,
                            moduleSize = input$ModuleSize_P6,
                            moduleSize_D2 = input$ModuleSize_D2_P6,
                            signed = input$signed_P6,
                            signed_D2 = input$signed_D2_P6,
                            cor_P3 = input$corNetwork_P6,
                            cor_D2_P3 = input$corNetwork_D2_P6,
                            cor_P4 = input$selectCorrelation_P6,
                            cor_D2_P4 = input$selectCorrelation_D2_P6,
                            loadExample = input$LoadExample2,
                            typeData_D1 = input$CountingT1,
                            typeData_D2 = input$OmicTable,
                            minimalCorrelation = input$minimalCorrelation)
    }
    
    
    file.copy("report.Rmd", tempReport(), overwrite = TRUE)
    
    
    
    # Set up parameters to pass to Rmd document
    if (input$LoadExample == "Yes"){
      params <- list(list_treatment = list_treatment,
                     list_input = list_input_P6,
                     exprDat = exprDat_report(),
                     sampleAnnot = sampleAnnot_Default(),
                     taxTable = taxTable_report())
    }else{
      if (input$LoadExample2 == "Yes"){
        params <- list(list_treatment = list_treatment,
                       list_input = list_input_P6,
                       exprDat = exprDat_report(),
                       sampleAnnot = sampleAnnot_Default(),
                       taxTable = taxTable_report(),
                       exprDatSec = exprDatSec_Default())
      }else{
        if (datasets[["type"]] == "S"){
          params <- list(list_treatment = list_treatment,
                         list_input = list_input_P6,
                         exprDat = exprDat_report(),
                         sampleAnnot = sampleAnnot())
        }else{
          if (datasets[["type"]] == "ST"){
            params <- list(list_treatment = list_treatment,
                           list_input = list_input_P6,
                           exprDat = exprDat_report(),
                           sampleAnnot = sampleAnnot(),
                           taxTable = taxTable_report())
          }else{
            if (datasets[["type"]] == "M"){
              params <- list(list_treatment = list_treatment,
                             list_input = list_input_P6,
                             exprDat = exprDat_report(),
                             sampleAnnot = sampleAnnot(),
                             exprDatSec = exprDatSec())
            }else{
              params <- list(list_treatment = list_treatment,
                             list_input = list_input_P6,
                             exprDat = exprDat_report(),
                             sampleAnnot = sampleAnnot(),
                             taxTable = taxTable_report(),
                             exprDatSec = exprDatSec())
            }
          }
        }
      }
    }
    
    
    incProgress(0.80, detail = "kniting report")
    paste(dirname(tempReport()), "report.html", sep = "")
    # Knit the document, passing in the `params` list, and eval it in a
    # child of the global environment (this isolates the code in the document
    # from the code in this app).
    rmarkdown::render(tempReport(), output_file = filename,
                      params = params,
                      envir = new.env(parent = globalenv()))
  }, message = 'Making report', value = 0)
})


output$downloadReport <- downloadHandler(
  
  # For PDF output, change this to "report.pdf"
  filename = paste(dirname(tempReport()), "/report.html", sep = ""),
  content = function(filename) {
    file.copy(paste(dirname(tempReport()), "/report.html", sep = ""), filename)
  }
)



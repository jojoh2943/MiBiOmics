#===============#
#### LIBRARY ####
#===============#

source("install_packages.R")


#### GENERAL FUNCTIONS ####

source("function.R")

#### GENERAL WGCNA CONFIGURATION ####
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#### JS AND CSS ####

css <- "
.nav li a.disabled {
  background-color: #aaa !important;
  color: #333 !important;
  cursor: not-allowed !important;
  border-color: #aaa !important;
}"


#==========#
#### UI ####
#==========#
# The UI generate the shiny web interface.

ui <- tagList(
  useShinyjs(),
  extendShinyjs(script = "function.js", functions = c("enableTab","disableTab")),
  inlineCSS(css),
  navbarPage(theme = shinytheme('yeti'),
             id = "navbar",
             "MiBiOmics",
             #### PAGE 1 ####
             # Page One allows the user to upload its files, filtrate, normalize, transform, remove outliers, view its table
             tabPanel("GENERAL_PARAMETERS", # The header tab Bar
                      id = "GENERAL_PARAMETERS",
                      icon = icon("folder"), #   adding icons --> go to https://fontawesome.com/icons?d=gallery
                      mainPanel(
                        fluidRow( 
                          column(6,
                                 tabsetPanel(type = "tabs", # A box containing several panels.
                                             tabPanel("Choose Your Analysis",
                                                      selectInput("TypeAnalysis",
                                                                  "Type of analysis: ",
                                                                  choices = c("Simple WGCNA" = "simple", "Multi-omics Analysis" = "multivariate", " " = "blank"),
                                                                  selected = "blank"),
                                                      conditionalPanel("input.TypeAnalysis == 'simple'",
                                                                       radioButtons("LoadExample", "Load the example dataset: ", choices = c("Yes"= "Yes", "No" = "No"), selected = "No"),
                                                                       conditionalPanel("input.LoadExample == 'No'",
                                                                                        radioButtons("CountingT", "Counting Table:", choices = c("Genes" = "Genes", "OTUs" = "OTUs"), selected = character(0)),
                                                                                        conditionalPanel("input.CountingT == 'OTUs'",
                                                                                                         checkboxInput("TaxonFile", "Additionnal Taxon File (CSV)", FALSE)))),
                                                      conditionalPanel("input.TypeAnalysis == 'multivariate'",
                                                                       radioButtons("LoadExample2", "Load the example dataset: ", choices = c("Yes"= "Yes", "No" = "No"), selected = "No"),
                                                                       conditionalPanel("input.LoadExample2 == 'No'",
                                                                                        radioButtons("CountingT1", "Counting Table:", choices = c("Genes" = "Genes", "OTUs" = "OTUs"), selected = character(0)),
                                                                                        conditionalPanel("input.CountingT1 == 'OTUs'",
                                                                                                         checkboxInput("TaxonFile1", "Additionnal Taxon File (CSV)", FALSE)),
                                                                                        radioButtons("OmicTable", "Second 'Omic' Table:", choices = c("Genes" = "Genes", "OTUs" = "OTUs", "Metabolites" = "Metabolites"), selected = character(0))
                                                                       )
                                                      )),
                                             tabPanel("Advanced Parameters",
                                                      radioButtons("Filtration", "Filtration:", choices = c("Yes" = "Yes", "No" = "No"), selected = "No"),
                                                      conditionalPanel("input.Filtration == 'Yes'",
                                                                       selectInput("TypeFiltration", "Type of filtration", choices = c("On prevalence" = "prevalence", "On minimum count" = "count")),
                                                                       conditionalPanel("input.TypeFiltration == 'prevalence'",
                                                                                        sliderInput("prevalence", "Prevalence Threshold: ", min = 0, max = 100, value = 0)),
                                                                       conditionalPanel("input.TypeFiltration == 'count'",
                                                                                        numericInput("count", "Minimum number of count across sample: ", min = 0, value = 2))),
                                                      hr(),
                                                      radioButtons("Normalisation", "Normalization", choices = c("Yes" = "Yes", "No" = "No"), selected = "No"),
                                                      conditionalPanel("input.Normalisation == 'Yes'",
                                                                       selectInput("TypeNormalisation", "Type of normalisation: ", choices = c("CSS" = "CSS", "TSS" = "TSS"))),
                                                      radioButtons("Transformation", "Transformation", choices = c("Yes" = "Yes", "No" = "No"), selected = "No"),
                                                      conditionalPanel("input.Transformation == 'Yes'",
                                                                       selectInput("TypeTransformation", "Type of transformation: ", choices = c("Log10" = "Log10", "Log2"= "Log2", "Hellinger"= "Hellinger", "Square" = "Square", "Square Root" = "sqrt", "ILR" = "ILR", "CLR" = "CLR"))),
                                                      hr(),
                                                      conditionalPanel("input.TypeAnalysis == 'multivariate'",
                                                                       radioButtons("Filtration1", "Filtration of the second table:", choices = c("Yes" = "Yes", "No" = "No"), selected = "No"),
                                                                       conditionalPanel("input.Filtration1 == 'Yes'",
                                                                                        selectInput("TypeFiltration1", "Type of filtration", choices = c("On prevalence" = "prevalence", "On minimum count" = "count")),
                                                                                        conditionalPanel("input.TypeFiltration1 == 'prevalence'",
                                                                                                         sliderInput("prevalence1", "Prevalence Threshold: ", min = 0, max = 100, value = 0)),
                                                                                        conditionalPanel("input.TypeFiltration1 == 'count'",
                                                                                                         numericInput("count1", "Minimum number of count across sample: ", min = 0, value = 2))),
                                                                       hr(),
                                                                       radioButtons("Normalisation1", "Normalization of the second table", choices = c("Yes" = "Yes", "No" = "No"), selected = "No"),
                                                                       conditionalPanel("input.Normalisation1 == 'Yes'",
                                                                                        selectInput("TypeNormalisation1", "Type of normalisation: ", choices = c("CSS" = "CSS", "TSS" = "TSS"))),
                                                                       radioButtons("Transformation1", "Transformation of the second table", choices = c("Yes" = "Yes", "No" = "No"), selected = "No"),
                                                                       conditionalPanel("input.Transformation1 == 'Yes'",
                                                                                        selectInput("TypeTransformation1", "Type of transformation: ", choices = c("Log10" = "Log10", "Log2"= "Log2", "Hellinger"= "Hellinger", "Square" = "Square", "Square Root" = "sqrt", "ILR" = "ILR", "CLR" = "CLR"))),
                                                                       hr()),
                                                      checkboxInput("batchEffect", "Batch Effect", value = FALSE),

                                                                       uiOutput("batchCol"),
                                                      uiOutput("selectSample")
                                             ) # end TabPanel advanced Parameter
                                 ) #End TabPanel where user select its options

                          ),
                          #### CHOOSING THE VIEWING OPTION
                          column(3,
                                 h3("Verify Your Tables"),
                                 uiOutput("ViewPanel")
                          ),
                          column(3,
                                 h3("Dimensions dataset:"),
                                 DT::dataTableOutput("ViewDim")
                                 )

                        ), #End first Row
                        fluidRow( #Begin Second Row --> Data visualization
                          h2("Information Visualization"),
                          DT::dataTableOutput("ViewTable")
                        ), #End second row

                        fluidRow( #Beginning third row

                          #### COUNTING TABLE UPLOAD
                          column(3,
                            h3("Upload Counting Table"),
                            conditionalPanel(("input.LoadExample == 'No' && input.LoadExample2 == 'No'"),
                                             fileInput("file1", "Choose CSV File (Counting Table): ",
                                                       multiple = FALSE,
                                                       accept = c("text/csv", "text/comma-separated-values,text/plain",".csv")),
                                             checkboxInput("header", "Header", TRUE),
                                             checkboxInput("rownames", "Rownames", TRUE),
                                             radioButtons("sep", "Separator", choices = c(Comma = ",",Semicolon = ";",Tab = "\t"), selected = ","),
                                             radioButtons("dec", "Decimal", choices = c(Comma = ",", Dot = "."), selected = ".")
                            )

                          ),

                          #### ANNOTATION TABLE UPLOAD
                          column(3,
                            h3("Upload Annotation File"),
                            conditionalPanel(("input.LoadExample == 'No' && input.LoadExample2 == 'No'"),
                                             fileInput("file2", "Choose CSV File (Annotation File): ",
                                                       multiple = FALSE,
                                                       accept = c("text/csv", "text/comma-separated-values,text/plain",".csv")),
                                             checkboxInput("header2", "Header", TRUE),
                                             checkboxInput("rownames2", "Rownames", TRUE),
                                             radioButtons("sep2", "Separator", choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), selected = ","),
                                             radioButtons("dec2", "Decimal", choices = c(Comma = ",", Dot = "."), selected = ".")
                            )
                          ),

                          #### TAXA TABLE UPLOAD
                          column(3,
                            h3("Upload Optional File"),
                            conditionalPanel(condition = "((input.CountingT == 'OTUs' | input.CountingT1 == 'OTUs') && (input.TaxonFile | input.TaxonFile1) && (input.LoadExample == 'No' | input.LoadExample2 == 'No'))",
                                             fileInput("file3", "Choose CSV File (Taxa Table): ",
                                                       multiple = FALSE,
                                                       accept = c("text/csv", "text/comma-separated-values,text/plain",".csv")),
                                             checkboxInput("header3", "Header", TRUE),
                                             checkboxInput("rownames3", "Rownames", TRUE),
                                             radioButtons("sep3", "Separator", choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), selected = ",")
                            )
                          ),

                          #### UPLOAD SECOND OMIC TABLE

                          column(3,
                            h3("Upload your second omic table"),
                            conditionalPanel("input.TypeAnalysis == 'multivariate' && input.LoadExample2 == 'No'",
                                             fileInput("file4", "Choose CSV File ('Omics' Counting Table): ",
                                                       multiple = FALSE,
                                                       accept = c("text/csv", "text/comma-separated-values,text/plain",".csv")),
                                             checkboxInput("header4", "Header", TRUE),
                                             checkboxInput("rownames4", "Rownames", TRUE),
                                             radioButtons("sep4", "Separator", choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), selected = ","),
                                             radioButtons("dec4", "Decimal", choices = c(Comma = ",", Dot = "."), selected = ".")
                            )
                          )

                      ) # End third row
                      )
             ),# End First Page

             #### PAGE 2 ####
             tabPanel("DATA_OVERVIEW",
                      id = "DATA_OVERVIEW",
                      icon = icon("signal", lib = "glyphicon"),
                      tabsetPanel(type = "tabs",
                                  id = "Page2_datasets",
                                  tabPanel("First Dataset", value = "Page2_dataset1",
                                           fluidRow( # Begin First row (preference selection)
                                             column(8,
                                                    tabsetPanel(type = "tabs", id = "Page2_dataset1",
                                                                tabPanel("Sample Dendrogramme", value = "Page2_dataset1_dendrogramme",
                                                                         plotOutput("Dendrogramme", width = 850, height = 650),
                                                                         selectInput("selectMethod",
                                                                                     label = "Choose a method for the dendrogramme: ",
                                                                                     choices = c("ward.D", "ward.D2", "single", "complete", "average", "median"),
                                                                                     selected = "average"),
                                                                         actionButton("showMethod", "How to choose the right method ?"),
                                                                         uiOutput("SelectVariable"),
                                                                         radioButtons("pdf_or_svg_page2_dataset1_dendrogramme",
                                                                                      "Choose an extension:",
                                                                                      choices = c("pdf", "svg")),
                                                                         downloadButton("Download_Page2_dataset1_dendrogramme", "Download")
                                                                ),
                                                                tabPanel("General Ordination", value = "Page2_dataset1_ordination",
                                                                         plotlyOutput("PCoA", width = 850, height = 650),
                                                                         selectInput("selectOrdination",
                                                                                     label = "Choose a type of ordination",
                                                                                     choices = c("PCA", "PCoA"),
                                                                                     selected = "PCA"),
                                                                         conditionalPanel("input.selectOrdination == 'PCoA'",
                                                                                          selectInput("selectDist",
                                                                                                      label = "Choose a distance for the PCoA: ",
                                                                                                      choices = c("bray", "euclidean", "manhattan", "jaccard", "binomial", "cao"),
                                                                                                      selected = "bray"),
                                                                                          actionButton("showDist", "How to choose the right distance ?")),
                                                                         uiOutput("SelectVariable3"),
                                                                         radioButtons("pdf_or_svg_page2_dataset1_ordination",
                                                                                      "Choose an extension:",
                                                                                      choices = c("pdf", "svg")),
                                                                         downloadButton("Download_Page2_dataset1_ordination", "Download")

                                                                ),
                                                                tabPanel("Dataset Specificity", value = "Page2_dataset1_rel_ab",
                                                                         plotOutput("Relative_Abundance", height = 650),
                                                                         uiOutput("SelectTaxo"),
                                                                         uiOutput("SelectVariable_rel_ab"),
                                                                         radioButtons("pdf_or_svg_page2_dataset1_rel_ab",
                                                                                      "Choose an extension:",
                                                                                      choices = c("pdf", "svg")),
                                                                         downloadButton("Download_Page2_dataset1_rel_ab", "Download")
                                                                ))
                                             )
                                  )),
                                  tabPanel("Second Dataset", value = "Page2_dataset2",
                                           fluidRow( # Begin First row (preference selection)
                                             column(8,
                                                    tabsetPanel(type = "tabs", id = "Page2_dataset2",
                                                                tabPanel("Sample Dendrogramme", value = "Page2_dataset2_dendrogramme",
                                                                         plotOutput("Dendrogramme_sec", width = 850, height = 650),
                                                                         selectInput("selectMethod_sec",
                                                                                     label = "Choose a method for the dendrogramme: ",
                                                                                     choices = c("ward.D", "ward.D2", "single", "complete", "average", "median"),
                                                                                     selected = "average"),
                                                                         actionButton("showMethod_sec", "How to choose the right method ?"),
                                                                         uiOutput("SelectVariable_sec"),
                                                                         radioButtons("pdf_or_svg_page2_dataset2_dendrogramme",
                                                                                      "Choose an extension:",
                                                                                      choices = c("pdf", "svg")),
                                                                         downloadButton("Download_Page2_dataset2_dendrogramme", "Download")
                                                                ),
                                                                tabPanel("General Ordination", value = "Page2_dataset2_ordination",

                                                                         plotlyOutput("PCoA_sec", width = 850, height = 650),
                                                                         selectInput("selectOrdinationSec",
                                                                                     label = "Choose a type of ordination",
                                                                                     choices = c("PCA", "PCoA"),
                                                                                     selected = "PCA"),
                                                                         conditionalPanel("input.selectOrdinationSec == 'PCoA'",
                                                                                          selectInput("selectDist_sec",
                                                                                                      label = "Choose a distance for the PCoA: ",
                                                                                                      choices = c("bray", "euclidean", "manhattan", "jaccard", "binomial", "cao"),
                                                                                                      selected = "bray"),
                                                                                          actionButton("showDist_sec", "How to choose the right distance ?")),

                                                                         uiOutput("SelectVariable3_sec"),
                                                                         radioButtons("pdf_or_svg_page2_dataset2_ordination",
                                                                                      "Choose an extension:",
                                                                                      choices = c("pdf", "svg")),
                                                                         downloadButton("Download_Page2_dataset2_ordination", "Download")

                                                                ),
                                                                tabPanel("Dataset Specificity", value = "Page2_dataset2_rel_ab",
                                                                         plotOutput("Relative_Abundance_sec", height = 650),
                                                                         uiOutput("SelectTaxo_sec"),
                                                                         uiOutput("SelectVariable_rel_ab_sec"),
                                                                         radioButtons("pdf_or_svg_page2_dataset2_rel_ab",
                                                                                      "Choose an extension:",
                                                                                      choices = c("pdf", "svg")),
                                                                         downloadButton("Download_Page2_dataset2_rel_ab", "Download")
                                                                )
                                                                )
                                             )
                                             ))


                      ), #End First Row
                      br()


  ), # End Second Page

  #### PAGE 3 ####
  tabPanel("NETWORK_INFERENCE", # Beginning third Page
           id = "NETWORK_INFERENCE",
           icon = icon("share-alt"),
           tabsetPanel(type = "tabs",
                       id = "NETWORK_INFERENCE_D1_D2",
                       tabPanel("First Dataset",
                                fluidRow( #Beginning first row
                                  column(4,
                                         h3("Choose Your Preferences: "),
                                         powers = c(c(1:10), seq(from = 12, to=20, by=2)),
                                         uiOutput("Power"),

                                         actionButton("showPower", "How to choose the right softPower ?"),

                                         hr(),

                                         uiOutput("ModuleSize"),

                                         hr(),
                                         radioButtons('signed', "Choose between a signed and a unsigned network: ",
                                                      choices = c("signed", "unsigned"),
                                                      selected = "signed"),

                                         actionButton("showSigned", "What is the difference between a signed and a unsigned Network ?"),

                                         hr(),

                                         radioButtons('corNetwork', "Choose the correlation option for the network: ",
                                                      choices= c('pearson', 'spearman', 'bicor'),
                                                      selected = 'pearson'),

                                         hr(),

                                         uiOutput("buttonPower"),

                                         uiOutput("buttonNet"),

                                         radioButtons("pdf_or_svg_p3",
                                                      "Choose an extension:",
                                                      choices = c("pdf", "svg")),

                                         downloadButton("Download_Network_Creation", "Download"),

                                         hr()

                                  ),


                                  column(6,
                                         h3("Scale Free Topoly against Soft Power: "),
                                         plotOutput("ScaleFree", height = 650)

                                  ),
                                  conditionalPanel("input.TypeAnalysis == 'simple'",
                                                   column(2,
                                                          h3("Click here to unlock the network exploration and the multivariate analysis
                                            (if you uploaded two counting tables)"),
                                                          actionButton("enableTabulations", "Explore your network")
                                                   )),
                                  conditionalPanel("input.TypeAnalysis == 'multivariate'",
                                                   column(2,
                                                          h3("Click here to unlock the network exploration and the multivariate analysis
                                            (if you uploaded two counting tables)"),
                                                          actionButton("goD2", "Check your second dataset")
                                                   ))


                                ), #End first Row
                                fluidRow( #Second Row

                                  column(6,
                                         h3("Classical Multidimensional Scaling plot to visualize the network"),
                                         uiOutput("multiScalePlotOutput")

                                  ),
                                  column(6,
                                         h3("Module Clustering: "),
                                         plotOutput("ModuleClustering")
                                  )
                                ), # End Second Row
                                fluidRow( #Beginning Third Row
                                  column(10,
                                         h3("Variables Dendrogramme with Module Colors: "),
                                         plotOutput("Dendrogramme_Colors"),
                                         downloadButton("downloadEigen", "Download Modules\' Eigenvalues")
                                  ),
                                  column(2,
                                         h3("Module Informations: "),
                                         DT::dataTableOutput("Module_Information"),
                                         downloadButton("downloadModule", "Download Modules")
                                  )
                                )
                                # , #End third Row
                                # h2("SparCC"),
                                # fluidRow( #Begining fourth row
                                #
                                #   column(10,
                                #          h3("SparCC Network"),
                                #          uiOutput("buttonSparcc"),
                                #          visNetworkOutput("sparcc.Network", height = 1000)
                                #   )
                                # ) #End Fourth Row
                       ),
                       tabPanel("Second Dataset", value = "NETWORK_INFERENCE_D2",
                                fluidRow( #Beginning first row
                                  column(4,
                                         h3("Choose Your Preferences: "),
                                         powers = c(c(1:10), seq(from = 12, to=20, by=2)),
                                         uiOutput("Power2"),

                                         actionButton("showPower2", "How to choose the right softPower ?"),

                                         hr(),

                                         uiOutput("ModuleSize2"),

                                         hr(),
                                         radioButtons('signed2', "Choose between a signed and a unsigned network: ",
                                                      choices = c("signed", "unsigned"),
                                                      selected = "signed"),

                                         actionButton("showSigned2", "What is the difference between a signed and a unsigned Network ?"),

                                         hr(),

                                         radioButtons('corNetwork2', "Choose the correlation option for the network: ",
                                                      choices= c('pearson', 'spearman', 'bicor'),
                                                      selected = 'pearson'),

                                         hr(),

                                         uiOutput("buttonPower2"),

                                         uiOutput("buttonNet2"),

                                         radioButtons("pdf_or_svg_p3_dataset2",
                                                      "Choose an extension:",
                                                      choices = c("pdf", "svg")),

                                         downloadButton("Download_Network_Creation_dataset2", "Download"),

                                         hr()

                                  ),
                                  column(6,
                                         h3("Scale Free Topoly against Soft Power: "),
                                         plotOutput("ScaleFree2", height = 650)

                                  ),
                                  column(2,
                                         h3("Click here to unlock the network exploration and the multivariate analysis
                                            (if you uploaded two counting tables)"),
                                         actionButton("enableTabulations2", "Explore your network")
                                  )

                                ), #End first Row
                                fluidRow( #Second Row

                                  column(6,
                                         h3("Classical Multidimensional Scaling plot to visualize the network"),
                                         uiOutput("multiScalePlotOutput2")

                                  ),
                                  column(6,
                                         h3("Module Clustering: "),
                                         plotOutput("ModuleClustering2")
                                  )
                                ), # End Second Row
                                fluidRow( #Beginning Third Row
                                  column(10,
                                         h3("Variables Dendrogramme with Module Colors: "),
                                         plotOutput("Dendrogramme_Colors2"),
                                         downloadButton("downloadEigenD2", "Download Modules\' Eigenvalues")
                                  ),
                                  column(2,
                                         h3("Module Informations: "),
                                         DT::dataTableOutput("Module_Information2"),
                                         downloadButton("downloadModuleD2", "Download Modules")
                                  )
                                )
                                # , #End third Row
                                # h2("SparCC"),
                                # fluidRow( #Begining fourth row
                                # 
                                #   column(10,
                                #          h3("SparCC Network"),
                                #          uiOutput("buttonSparcc2"),
                                #          visNetworkOutput("sparcc.Network2", height = 1000)
                                #   )
                                # ) # End Fourth Row
                       )
                       )
           ), # End third Page

 ### PAGE 4 ###
  # tabPanel("SPARCC", #Beginning page 4
  #          id = "SPARCC",
  #          icon = icon("spinner"),
  #          fluidRow( #First Row
  #            column(9,
  #                   h3("SparCC Network"),
  #                   visNetworkOutput("sparcc.Network")
  #            )
             # ,
             # column(3,
             #        h3("Module information:"),
             #        DT::dataTableOutput("sparcc_table")
             # )
           # )
           # ,# End First Row
           # fluidRow( #beginning second row
           #   column(3,
           #          h3("Choose your preferences:"),
           #          uiOutput("SelectModule_sparcc"),
           #          uiOutput("SelectTaxo_sparcc"))
           # ), #End second row
           # fluidRow(
           #   column(6,
           #          h3("Relative Abundance SparCC Network Modules"),
           #          #verbatimTextOutput("Relative_Abundance_Module_sparcc")
           #          plotlyOutput("Relative_Abundance_Module_sparcc", height = 650)
           #
           #   ),
           #   column(6,
           #          h3("Correlation between the modules and the external variables"),
           #          #verbatimTextOutput("Relative_Abundance_Module_sparcc")
           #          iheatmaprOutput("HEATMAP_SPARCC", height = 650)
           #
           #   )
           # ) #Beginning third row

  # ), # End Page 5

 #### PAGE 4 ####
 tabPanel("NETWORK_EXPLORATION", # Beginning fourth page
          id = "NETWORK_EXPLORATION",
          icon = icon("share-alt-square"),
          tabsetPanel(type = "tabs",
                      tabPanel("First Dataset",
                               fluidRow( # first row
                                 column(3,
                                        h3("Choose Your Preferences"),
                                        p("In this section you will be able to explore the genes or OTUs of specific modules."),
                                        p("The selection inputs under each plots allow you to select your preferences"),
                                        p("First, take a look at the heatmap on your right and try to find interesting correlation. For the
                                          modules and variables of your choice, explore their specificity in the network to better understand what
                                          links these entity together."),
                                        hr(),
                                        radioButtons("selectCorrelation", label = "Choose a correlation method to link modules to external traits:",
                                                     choices = list("Spearman" = "spearman", "Pearson" = "pearson", "Kendall" = "kendall"),
                                                     selected = "spearman"),
                                        uiOutput("SelectModule"),
                                        hr(),
                                        
                                        radioButtons("pdf_or_svg_p4",
                                                     "Choose an extension:",
                                                     choices = c("pdf", "svg")),
                                        sliderInput("widthPDF", "Figures Width:", min = 10, max = 100, value = 30),
                                        sliderInput("heightPDF", "Figures Height: ", min = 10, max = 100, value = 30),
                                        downloadButton("Download_Network_Exploration", "Download")
                                        
                                        ),
                                 column(9,
                                        title = "Module's Correlation to External Traits",
                                        plotOutput("Corr_External_Trait", height = 750)
                                 )
                                 
                                 ), # end first row
                               br(),
                               br(),
                               fluidRow( # Beginning Second Row
                                 column(6,
                                        h3("Contribution of the Samples to the Modules"),
                                        plotOutput("Sample_Contribution", height = 650),
                                        uiOutput("SelectVariable_barplot")
                                 ),
                                 column(6,
                                        h3("Variable's Correlation to External Trait Against Module Membership"),
                                        plotOutput("Module_Membership", height = 650),
                                        uiOutput("SelectVariable1")
                                 )
                               ), # End Second Row
                               br(),
                               hr(),
                               br(),
                               fluidRow( # Beginning Third Row
                                 column(9,
                                        title = "Module Specificity",
                                        plotOutput("Relative_Abundance_Module", height = 1000)
                                 ),
                                 column(3,
                                        uiOutput("SelectTaxo1"),
                                        uiOutput("SelectVariable_RelAb")
                                 )
                               ),
                               br(),
                               fluidRow(
                                 h2("PLS and VIP Scores"),
                                 column(3,
                                        uiOutput("sampleAnnotSelection"),
                                        uiOutput("colorModule")),
                                 column(3,
                                        plotOutput("ncomp"),
                                        numericInput("ncomponent", "Number of Components:",
                                                     1, min = 1, max = 10,
                                                     value = 6),
                                        
                                        
                                        downloadButton("PLS_VIP", "Download")),
                                 column(6,
                                        h3("PLS"),
                                        plotOutput("PLS"),
                                        h5("CVS:"),
                                        textOutput("CVS"))
                               ),
                               fluidRow(
                                 column(9,
                                        h2("Hive Plot"),
                                        uiOutput("taxAnnot"),
                                        plotOutput("edge_node", height = 1000),
                                        radioButtons("pdf_or_svg_hive_p4",
                                                     "Choose an extension:",
                                                     choices = c("pdf", "svg")),
                                        sliderInput("HIVEwidthPDF", "Figures Width:", min = 10, max = 100, value = 30),
                                        sliderInput("HIVEheightPDF", "Figures Height: ", min = 10, max = 100, value = 30),
                                        downloadButton("hivePlot_D1", "Download")
                                        
                                 )
                               )
                               
                               # End Third Row
          ),
          tabPanel("Second Dataset",
                   fluidRow( # first row
                     column(3,
                            h3("Choose Your Preferences"),
                            p("In this section you will be able to explore the genes or OTUs of specific modules."),
                            p("The selection inputs under each plots allow you to select your preferences"),
                            p("First, take a look at the heatmap on your right and try to find interesting correlation. For the
                              modules and variables of your choice, explore their specificity in the network to better understand what
                              links these entity together."),
                            hr(),
                            radioButtons("selectCorrelationSec", label = "Choose a correlation method to link modules to external traits:",
                                         choices = list("Spearman" = "spearman", "Pearson" = "pearson", "Kendall" = "kendall"),
                                         selected = "spearman"),
                            uiOutput("SelectModuleSec"),
                            hr(),
                            radioButtons("pdf_or_svg_p4_dataset2",
                                         "Choose an extension:",
                                         choices = c("pdf", "svg")),
                            sliderInput("widthPDF_D2", "Figures Width: ", min = 10, max = 100, value = 30),
                            sliderInput("heightPDF_D2", "Figures Height: ", min = 10, max = 100, value = 30),
                            
                            downloadButton("Download_Network_Exploration_dataset_2", "Download")
                            
                            )
                     ,
                     column(9,
                            title = "Module's Correlation to External Traits",
                            plotOutput("Corr_External_TraitSec", height = 1000, width = 1000)
                     )
                     
                     ), # end first row
                   br(),
                   br(),
                   fluidRow( # Beginning Second Row
                     column(6,
                            h3("Contribution of the Samples to the Modules"),
                            plotOutput("Sample_ContributionSec", height = 650),
                            uiOutput("SelectVariable_barplotSec")
                     ),
                     column(6,
                            h3("Variable's Correlation to External Trait Against Module Membership"),
                            plotOutput("Module_MembershipSec", height = 650),
                            uiOutput("SelectVariable1Sec")
                     )
                   ), # End Second Row
                   br(),
                   hr(),
                   br(),
                   fluidRow( # Beginning Third Row
                     column(6,
                            title = "Module Specificity",
                            plotOutput("Relative_Abundance_ModuleSec", height = 650)
                     ),
                     column(3,
                            uiOutput("SelectTaxo1Sec"),
                            uiOutput("SelectVariable_RelAbSec")
                     )
                   ),
                   br(),
                   fluidRow(
                     h2("PLS and VIP Scores"),
                     column(3,
                            uiOutput("sampleAnnotSelection_D2"),
                            uiOutput("colorModule_D2")),
                     column(3,
                            plotOutput("ncomp_D2"),
                            numericInput("ncomponent_D2", "Number of Components:",
                                         1, min = 1, max = 10,
                                         value = 6),
                            
                            downloadButton("PLS_VIP_D2", "Download")),
                     column(6,
                            h3("PLS"),
                            plotOutput("PLS_D2"),
                            h5("CVS:"),
                            textOutput("CVS_D2"))
                     
                   ),
                   fluidRow(
                     column(9,
                            h2("Hive Plot"),
                            uiOutput("taxAnnotSec"),
                            plotOutput("edge_node_D2", height = 1000),
                            radioButtons("pdf_or_svg_hiveD2_p4",
                                         "Choose an extension:",
                                         choices = c("pdf", "svg")),
                            sliderInput("HIVED2widthPDF", "Figures Width:", min = 10, max = 100, value = 30),
                            sliderInput("HIVED2heightPDF", "Figures Height: ", min = 10, max = 100, value = 30),
                            downloadButton("hivePlot_D2", "Download")
                            
                     )
                     
                   )
 )
          )
                      ), # End fourth page
 

  #### PAGE 5 ####
  tabPanel("MULTI-OMICS_ANALYSIS", #Beginning page 5
           id = "MULTI-OMICS_ANALYSIS",
           icon = icon("sitemap"),
           fluidRow( #First Row
             column(3,
               h3("Select your preferences: "),
               h5("Choose an option for the first dataset:"),
               uiOutput("SelectModule1"),

               hr(),

               h5("Choose an option for the first dataset:"),
               uiOutput("SelectModule2"),

               hr(),

               uiOutput("SelectVariable2"),

               hr(),

               uiOutput("SelectTaxonomy"),

               hr(),

               checkboxInput("ShowDrivers", "Show the main drivers of the co-inertia", FALSE),
               downloadButton("Download_Coinertia_drivers", "Download drivers"),

               radioButtons("pdf_or_svg_p5",
                            "Choose an extension:",
                            choices = c("pdf", "svg")),

               downloadButton("Download_Multivariate_Analysis", "Download")
             ),
             column(9,
                    tabsetPanel(type = "tabs",
                                tabPanel("Co-Inertia Analysis",
                                         plotOutput("coinertia", height = 750, width = 750),
                                         h4("RV-score: "),
                                         verbatimTextOutput("RV")
                                ),
                                tabPanel("Co-inertia Information",
                                         plotOutput("coinertia_bivariate", height = 750),
                                         h4("RV-score: "),
                                         verbatimTextOutput("RV2")
                                         ),
                                tabPanel("Procrustes Analysis",
                                         plotOutput("PA", height = 750)
                                         )
                    )
             )
           ), # End First Row
           fluidRow(12, # Second Row
                    h3("General Heatmap representing the correlation between the modules eigengenes of each network."),
                    iheatmaprOutput("HEATMAP_MEs", height = 1000)),
           fluidRow( # Third Row
             column(6,
                    h3("Choose an option for the first dataset"),
                    uiOutput("SelectModule11")),
                    conditionalPanel("input.selectModule11_p5 == 'Individual_Variables'",
                              h4("Select some variables to show in the heatmap: "),
                              uiOutput("selectVariablesExprDat")),
             column(6,
                    h3("Choose an option for the second dataset"),
                    uiOutput("SelectModule22"),
                    conditionalPanel("input.selectModule22_p5 == 'Individual_Variables'",
                                     h4("Select some variables to show in the heatmap: "),
                                     uiOutput("selectVariables"))

             )
           ), # End Third Row
           fluidRow( #Fourth Row
             column(12,
                    h3("Bipartite_network representing the associations between modules of each networks"),
                    plotOutput("bipartite_network", height = 1000),
                    h3("Heatmap representing the correlation between the module of the network and the selected variables of the second dataframe uploaded"),
                    iheatmaprOutput("HEATMAP", height = 1000)
             )
           ) #End Fourth Row
           ), # End Page 5

 #### PAGE 6 ####
 tabPanel("REPORT_GENERATION",
          fluidRow(
            column(12,
                   h3("Generate the report of all the results:"),
                   h4("Section 1: Data import and treatment"),
                   textAreaInput("commentSection1", "Add a comment about the data import and treament (First Section)", width = "700px", height = "200px"),
                   hr(),
                   h4("Section 2: Data overview"),
                   selectInput("selectMethod_P6",
                               label = "Choose a method for the sample dendrogramme: ",
                               choices = c("ward.D", "ward.D2", "single", "complete", "average", "median"),
                               selected = "average"),
                   conditionalPanel("input.TypeAnalysis == 'multivariate'",
                                    selectInput("selectMethod_D2_P6",
                                                label = "Choose a method for the sample dendrogramme (second Dataset): ",
                                                choices = c("ward.D", "ward.D2", "single", "complete", "average", "median"),
                                                selected = "average")),
                   actionButton("showMethod_P6", "How to choose the right method ?"),
                   selectInput("selectOrdination_P6",
                               label = "Choose a type of ordination",
                               choices = c("PCA", "PCoA"),
                               selected = "PCA"),
                   conditionalPanel("input.selectOrdination_P6 == 'PCoA'",
                                    selectInput("selectDist_P6",
                                                label = "Choose a distance for the PCoA: ",
                                                choices = c("bray", "euclidean", "manhattan", "jaccard", "binomial", "cao"),
                                                selected = "bray")),
                   conditionalPanel("input.TypeAnalysis == 'multivariate'",
                                    selectInput("selectOrdination_D2_P6",
                                                label = "Choose a type of ordination",
                                                choices = c("PCA", "PCoA"),
                                                selected = "PCA"),
                                    conditionalPanel("input.selectOrdination_D2_P6 == 'PCoA'",
                                                     selectInput("selectDist_D2_P6",
                                                                 label = "Choose a distance for the PCoA: ",
                                                                 choices = c("bray", "euclidean", "manhattan", "jaccard", "binomial", "cao"),
                                                                 selected = "bray"))),
                   actionButton("showDist_P6", "How to choose the right distance ?"),
                   textAreaInput("commentSection2", "Add a comment about the data overview (Second Section)", width = "700px", height = "200px"),
                   hr(),

                   h4("Section 3: Network Creation"),
                   column(6,
                          h4("FIRST DATASET:"),
                          h6("Scale Free Topoly against Soft Power (First dataset): "),
                          plotOutput("ScaleFree_P6", height = 650),
                          uiOutput("Power_P6"),
                          actionButton("showPower_P6", "How to choose the right softPower ?"),
                          uiOutput("ModuleSize_P6"),
                          radioButtons('signed_P6', "Choose between a signed and a unsigned network (First dataset): ",
                                       choices = c("signed", "unsigned"),
                                       selected = "signed"),
                          actionButton("showSigned_P6", "What is the difference between a signed and a unsigned Network ?"),
                          radioButtons('corNetwork_P6', "Choose the correlation option for the network (First dataset): ",
                                       choices= c('pearson', 'spearman', 'bicor', 'sparCC'),
                                       selected = 'pearson')),
                   conditionalPanel("input.TypeAnalysis == 'multivariate'",
                     column(6,
                          h4("SECOND DATASET:"),
                          h6("Scale Free Topoly against Soft Power (Second dataset): "),
                          plotOutput("ScaleFree_D2_P6", height = 650),
                          uiOutput("Power_D2_P6"),
                          actionButton("showPower_D2_P6", "How to choose the right softPower ?"),
                          uiOutput("ModuleSize_D2_P6"),
                          radioButtons('signed_D2_P6', "Choose between a signed and a unsigned network (Second dataset): ",
                                       choices = c("signed", "unsigned"),
                                       selected = "signed"),
                          actionButton("showSigned_D2_P6", "What is the difference between a signed and a unsigned Network ?"),
                          radioButtons('corNetwork_D2_P6', "Choose the correlation option for the network (Second dataset): ",
                                       choices= c('pearson', 'spearman', 'bicor', 'sparCC'),
                                       selected = 'pearson'))
                     ),
                   textAreaInput("commentSection3", "Add a comment about the network inference (Third Section)", width = "700px", height = "200px"),
                   hr(),
                   h4("Section 4: Network Exploration"),
                   column(6,
                          h4("FIRST DATASET:"),
                          radioButtons("selectCorrelation_P6", label = "Choose a correlation method for the first dataset \nto link modules to external traits:",
                                       choices = list("Spearman" = "spearman", "Pearson" = "pearson", "Kendall" = "kendall"),
                                       selected = "spearman")),
                   conditionalPanel("input.TypeAnalysis == 'multivariate'",
                                    column(6,
                                           h4("SECOND DATASET"),
                                           radioButtons("selectCorrelation_D2_P6", label = "Choose a correlation method for the first dataset \nto link modules to external traits:",
                                                        choices = list("Spearman" = "spearman", "Pearson" = "pearson", "Kendall" = "kendall"),
                                                        selected = "spearman"))),
                   textAreaInput("commentSection4", "Add a comment about the network exploration (Fourth Section)", width = "700px", height = "200px"),
                   hr(),
                   conditionalPanel("input.TypeAnalysis == 'multivariate'",
                                    h4("Section 5: Multi-omics analysis:"),
                                    numericInput("minimalCorrelation", "Choose a minimum correlation input to determine which figures will be printed in the network: ", min = 0.50, max = 1.00, step = 0.01, value = 0.50),
                                    textAreaInput("commentSection5", "Add a comment about the multi-omics analysis (Fith Section)", width = "700px", height = "200px")),

                   actionButton("produceReport", "Produce a Report"),
                   downloadButton("downloadReport", "Report Generation"))
          )),

  #### HELP PAGE ####
  tabPanel("HELP", # Beginning Help Page
           fluidRow(
             column(12,
               h3("What is the purpose of Network Explorer ?"),
               style = "background-color: #fcd688; box-shadow: 5px 10px;",
               p("MiBiOmics uses several R libraries in a user friendly way in order to explore and to understand the correlation patterns
                                 in one or several omics datasets. Depending on the analysis you seek, you will be able to study a single dataset through
                                 WGCNA correlation network, ordination techniques and dendrogrammes to find what drives the overall variability of your
                                 data but you will also be able to compare two omic datasets in what we call 'multi-omics analysis' using basic
                                 correlation, network analysis, co-inertia and procrustes analysis."),
               br(),
               p("We propose here a methodology for a system biology approach where it is possible to see, compare and understand the
                                 different aspect of the same biological object. More and more, we see in experiments, several types of omics data describing
                                 the same individuals. However, it remains really difficult to compare them, even if they are part of the same system. We know
                                 today that genes, RNAs, proteins, OTUs, metabolites, are not closed systems. They interact with each other. The purpose of
                                 our method is to find the associations between these different entities.")
             )
           ),
           hr(),
           fluidRow(
             column(4,
                    style = "background-color:   #f17e70 ; box-shadow: 5px 10px;",
                    h3("First Step: Upload your datasets ! Filtrate ! Normalize ! Transform !"),
                    h4("If you want to perform a simple WGCNA:"),
                    p("You only have to upload two datasets. The first one is your counting table. The counting table contains the
                                 name of your genes/OTUs in the columns and the names of your samples in the rows. Each cells corresponds to the
                                 number of the correponding gene/OTU found in this sample. (see Table 1)"),
                    p("The counting table must be accompagnied by a Annotation table containing the exact same sample names as rows
                                 and several descriptive variables as columns (it can be day of experiments, special treatment, collection site ...). (see Table 2)"),
                    p("If your working with OTUs, you can also upload a Taxa Table describing, for each OTUs, the corresponding Phylum,
                                 Order, Class, ..., Species but selecting the option 'upload a taxa table'. (see Table 3)"),
                    p("In advance parameters, you can filtrate, normalize, transform your data. You can also remove sample outliers."),
                    h4("If you want to perform a multi-omics analysis:"),
                    p("MiBiOmics allows to perform analysis on two omics layers to compare, confront and associate different variables describind the same individuals.
                      To perform a multi-omics analysis, select the Multi-Omics Analysis fiels in Type of analysis. This time you will have to upload at least 2 counting
                      tables (the first one contains the variables of your first omics layer and the second one contains the variables of your second omics layer). You will
                      also need, as before, a annotation table. If one of your Omics layer contains OTUs, you need to upload this counting table first to be able to upload its
                      associated taxonomic annotation. As before, you will also be able to filtrate, normalize and transform both counting tables.")
             ),
             column(4,
               DT::dataTableOutput("exampleTable")
             ),
             column(3,
               DT::dataTableOutput("exampleAnnot")
             ),
             column(8,
                    DT::dataTableOutput("exampleTaxa")
             )
           ),

           hr(),
           fluidRow(
             column(5,
                    style = "background-color: #86c5d6; box-shadow: 5px 10px;",
                    h3("Explore your dataset !"),
                    p("When you are ready, you uploaded your datasets, you treated your data with filtration and/or normalization
                                 and/or transformation, you will be able to explore your dataset(s) variability in the second page of MiBiOmics. In this section, a PCA or a PCoA (depending on if you are working with genes
                                 or with OTUs) will allow you to see the main axis of variability. You can change the color of your samples
                                 according to the external variables uploaded with the annotation table to see what might play a part in this
                                 variability."),
                    p("A cluster dendrogramme is also plotted to help you visualize which samples are closely related."),
                    p("If you uploaded a OTUs Counting Table with the corresponding Taxa Table you will be able to see a figure
                                 representing the relative taxa abundance of each sample at different taxonomic level. "),
                    h4("If you perform a multi-omics analysis"),
                    p("Use the tabset called 'second dataset' to explore the variability of your second omics layer."),
                    br()
             ),
             column(5,
                    offset = 2,
                    style = "background-color: #a987d6; box-shadow: 5px 10px;",
                    h3("Create a correlation Network with WGCNA !"),
                    p("WGCNA is a R package which uses a system biology approach to find patterns of correlation in omics datasets (Peter Langfelder
                                 and Steve Horvath, 2008). Based on correlation networks, the method defines modules as groups of highly interconnected
                                 genes and reduce the dimensionality of the data by helping the user concentrate on the module of interest. A motivation
                                 for the implementation of this application was to give the ability to create the network in a user-friendly way and to
                                 understand better the influence of the principal parameters, like the soft power or the minimum module size, on the
                                 structure of the Network."),
                    p("Once your in the
                                 Network Inference section, you will be able to change the soft power and minimum size module to see how the structure of your
                                 correlation network is affected. This step migth take a while to be completed if the counting table uploaded is too large ( > 5000).
                                 Unfortunately, it is not possible to perform correlation network with datasets larger than 10 000 and we are currently working on
                                 this issue."),
                    h4("If you perform a multi-omics analysis"),
                    p("Use the tabset called 'second dataset' to perform the network inference on your second omics layer. /!\ The next section, 'Network exploration' can
                      only be used once the network inference step has been performed on both datasets.")
                    )
           ),
           hr(),
           fluidRow(
             column(12,
                    style = "background-color:   #ff9900; box-shadow: 5px 10px;",
                    h3("Explore your Network !"),
                    p("The Network exploration page allows you to characterise the modules, to understand what are the common properties
                               of the genes/OTUs belonging to a same group. In this section you will be able to see which samples contribute the
                               most to the formation of a module, how modules correlate to external traits and the relationship between the module
                               membership of a gene/OTU and its correlation to external traits. In this page you can explore the inside of the
                               network, part by part, and associate a group of highly interconnected genes/OTUs to a particular aspect of the
                               experiment."),
                    p("If you uploaded an OTUs counting table, you will be able to explore the relative taxon abundance of each module."),
                    h4("If you perform a multi-omics analysis"),
                    p("Use the tabset called 'second dataset' to explore the modules of your second network."),
                    p("This section is important to select the modules of interest in each omics layers. To compare modules of the two different networks, they need to
                      to be associated to the same external traits. For example the module yellow of your first dataset and the module blue of your second dataset are
                      strongly associated to the sampleSite external trait. They might be used and associated in the multi-omics analysis tab of MiBiOmics.")
             )

           ),
           hr(),
           fluidRow(
             column(12,
                    style = "background-color: #fcd688; box-shadow: 5px 10px;",
                    h3("Perform a multi-omics analysis (only if the option multi-omics analysis was selected in the first tab) !"),
                    h4("Ordination technique: Co-Inertia Analysis"),
                    p("The co-inertia analysis is performed with the omicade4 R package used to perform multiple co-inertia analysis and created by Meng C, et al. in 2013.
                      This package is based on the work realized by Stephane Dray, Daniel Chessel and Jean Thioulouse in 2003. With the co-inertia they created a methodoly
                      to associate two datasets and to represent the co-variance between both datasets. In MiBiOmics, the co-inertia is represented in two different ways.
                      In the first figure (on the top in the multi-omics analysis tab), the samples are placed on the plot according to their values in the co-inertia first
                      and second axis (read the article Dray et al., 2003 for more details). The closest they are to the maximal correlation line in red, the more they are
                      correlated according to their expreesions/concentrations/abundances values in both datasets. The second figure represents the co-variance between the
                      values of each sample in both datasets. Each sample is represented twice in the plot: the filled circle represents the first datasets and more specifically
                      the position of the sample in the first ordination space, the empty square represents the second datasets and the position of the same sample in the second
                      ordination space. Both representation of the samples are linked with a line and the length of the line indicates how much the values of the same sample covary
                      from one dataset to the over. The longer is the line the less the values of the same sample covary from one dataset to the other. The variables that are the most
                      associated with the covariance of both datasets are also plotted on the co-inertia. They can belong to both datasets. The RV score indicates the quality
                      of the co-inertia. It values varies between 0 and 1 and the closest it is to 1, the best is the covariation between both datasets."),
                    p("The co-inertia information tabset shows the length of this line for each samples represented with a dendrogramme, and two other figures. This tabsets might help
                      to isolate the samples with a low covariation and see if it can be associated with a external clinical traits."),
                    h4("Ordination Technique: Procrustes Analysis"),
                    p("Procrustes analysis is realized in MiBiOmics thanks to the vegan R package created by Jari Oksanen. The procrustes method was originally created
                      by J. C. Gower in 1975 and uses the 3D configuration of the datasets ordinations, reorientate and re-scale them to allow the comparison of these differents
                      ordination in multidimensional space. The procrustes analysis plot can be found in the procrustes analysis tabset and can be interpreted as the co-inertia
                      analysis plot."),
                    h4("Correlation Network analysis: using WGCNA modules to compare sub-parts of two different networks."),
                    p("First, the first heatmap represents the correlation between each modules of both networks. The more two modules are correlated the more they might be
                      associated to each other. This heatmap needs to be used as a guide to select the module of both network (and both omics layers) to each other. Only one
                      module of each network can be selected to pursue the analysis."),
                    p("Once a module from each network (each omics dataset) is selected, the pairwise correlation between each variables of each modules is performed and represented
                      in two different manners: a bi-partite networks and a correlation heatmap. In the bi-partite networks, red nodes represents the variables of the first dataset
                      and blue nodes the variables of the second datasets. A line between two nodes indicates an association. The following heatmap represents the correlation values
                      between the variables of the first selected module (first omics layer) and the variables of the second selected module (second omic layer)")
                    ),
             hr()

                    )
           ), # End Help Page
 #### ABOUT ####
 tabPanel("ABOUT",
          fluidRow(
            column(6, offset = 1,
                   h2("Citations"),
                   p("MiBiOmics v0.0.1 was developed with Shiny, a RStudio package for web application development. The methodology available on the application relies
                     on the following R packages: "),
                   a(href="https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/", "WGCNA"),
                   p(""),
                   a(href="https://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf", "vegan"),
                   p(""),
                   a(href="https://www.bioconductor.org/packages/release/bioc/vignettes/omicade4/inst/doc/omicade4.pdf", "omicade4"),
                   p("For WGCNA, we strongly advice to read their publication and look at their tutorial to choose the parameter in the network inference section:"),
                   a(href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559", "Peter Langfelder and Steve Horvath, 2008"),
                   p("Procrustes and co-inertia analysis were performed with vegan and omicade4 functions respectively. Please refer to the original publication
                     to aknowledge the methodology:"),
                   a(href="https://link.springer.com/content/pdf/10.1007%2FBF02291478.pdf", "J. C. Gower, 1975 (procrustes analysis)"),
                   p(""),
                   a(href="http://pbil.univ-lyon1.fr/members/dray/files/articles/dray2003c.pdf", "S. Dray et al., 2003 (co-inertia analysis)")
                   )

          )
 ) # end about page

)
)

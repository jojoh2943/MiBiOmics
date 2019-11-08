# Define server logic ----
server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  options(warn=1)
  
  #### SERVER FUNCTIONS ####
  source("server_function/SERVER_Data_Upload.R", local=TRUE)
  source("server_function/SERVER_Data_Exploration.R", local=TRUE)
  source("server_function/SERVER_Network_Inference.R", local=TRUE)
  source("server_function/SERVER_Network_Exploration.R", local=TRUE)
  source("server_function/SERVER_Multi_Omics_Analysis.R", local=TRUE)
  source("server_function/SERVER_Report_Help.R", local=TRUE)

  #### GENERAL CODE ####

  #### HEADER ####  
  js$disableTab("NETWORK_EXPLORATION")
  js$disableTab("MULTI-OMICS_ANALYSIS")
  
  observeEvent(input$enableTabulations, {
    # enable NETWORK_EXPLORATION when clicking the button
    js$enableTab("NETWORK_EXPLORATION")
    #updateTabsetPanel(session, "NETWORK_INFERENCE_D1_D2", "NETWORK_INFERENCE_D2")
    # switch to NETWORK_EXPLORATION
    updateTabsetPanel(session, "navbar", "NETWORK_EXPLORATION")
  })
  
  observeEvent(input$goD2, {
    updateTabsetPanel(session, "NETWORK_INFERENCE_D1_D2", "NETWORK_INFERENCE_D2")
  })
  
  observeEvent(input$goD3, {
    updateTabsetPanel(session, "NETWORK_INFERENCE_D1_D2", "NETWORK_INFERENCE_D3")
  })
  
  observeEvent(input$enableTabulations_D3, {
    # enable NETWORK_EXPLORATION when clicking the button
    js$enableTab("NETWORK_EXPLORATION")

      js$enableTab("MULTI-OMICS_ANALYSIS")
  
    # switch to NETWORK_EXPLORATION
    updateTabsetPanel(session, "navbar", "NETWORK_EXPLORATION")
  })
  
  observeEvent(input$enableTabulations_D2, {
    # enable NETWORK_EXPLORATION when clicking the button
    js$enableTab("NETWORK_EXPLORATION")
    
    js$enableTab("MULTI-OMICS_ANALYSIS")
    
    # switch to NETWORK_EXPLORATION
    updateTabsetPanel(session, "navbar", "NETWORK_EXPLORATION")
  })
  
 
}






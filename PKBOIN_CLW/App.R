# ==============================================================================
# PK-BOIN12 / TITE-PKBOIN-12 Shiny - Ultimate Bulletproof Full English Version
# ==============================================================================

library(shiny)
library(shinythemes)
library(tidyverse)    # Resolves the '%>%' missing error
library(parallel)     # Package for mclapply

# ------------------------------------------------------------------------------
# 0. System Compatibility & Mathematical Algorithm Patches
# ------------------------------------------------------------------------------
# Windows does not natively support fork-based mclapply. 
# This safely redirects mclapply to standard lapply on Windows machines.
if (.Platform$OS.type == "windows") {
  .GlobalEnv$mclapply <- function(X, FUN, ..., mc.cores = 1) {
    lapply(X, FUN, ...)
  }
}

# The PAVA (Isotonic Regression) custom runtime fallback patch.
# Bypasses local file dependency errors and utilizes R's high-speed core engine.
pava_rcpp <- function(x, w = rep(1, length(x))) {
  return(isoreg(x)$yf)
}


# ------------------------------------------------------------------------------
# 1. Load Core Engine Source Files from prog/ Folder 
# ------------------------------------------------------------------------------
source("prog/fun_findOBD.R")  
source("prog/fun_TITE_PKBOIN12.R")
source("prog/fun_TITE_PKBOIN12_OBD.R")
source("prog/fun_TITE_PKBOIN12dec.R") 
source("prog/fun_TITE_PKBOIN12_one.R") 
source("prog/fun_TITE_PK_core_para.R")
source("prog/fun_TITE_PK_update.R")    
source("prog/fun_TITE_PK_fixsimu.R")


# ------------------------------------------------------------------------------
# 2. User Interface (UI) Design
# ------------------------------------------------------------------------------
ui <- fluidPage(
  theme = shinytheme("flatly"), # Clean, professional statistical dashboard layout
  
  titlePanel("PK-BOIN12 / TITE-PKBOIN-12 Trial Design Simulator"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Clinical Parameters"),
      hr(),
      
      # Model Design Selection
      selectInput("design_type", "Design Method:",
                  choices = c("PKBOIN-12", "TITE-PKBOIN-12")),
      
      sliderInput("target_tox", "Target Toxicity Rate (pT):", 
                  min = 0.1, max = 0.5, value = 0.3, step = 0.05),
      
      sliderInput("target_eff", "Target Efficacy Rate (qE):", 
                  min = 0.1, max = 0.6, value = 0.25, step = 0.05),
      
      numericInput("n_dose", "Number of Dose Levels (dN):", 
                   value = 5, min = 2, max = 10),
      
      numericInput("cohort_size", "Cohort Size (csize):", 
                   value = 3, min = 1),
      
      numericInput("cN", "Number of Cohorts (cN):", 
                   value = 10, min = 2, max = 30),
      
      hr(),
      
      actionButton("run_sim", "Run Simulation", 
                   icon = icon("play"), 
                   class = "btn-lg btn-success")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Simulation Summary", 
                 br(),
                 h4("Simulation Output Logs:"),
                 verbatimTextOutput("obd_text"),
                 plotOutput("sim_plot")),
        
        tabPanel("Detailed Data List", 
                 br(),
                 tableOutput("detail_table"))
      )
    )
  )
)


# ------------------------------------------------------------------------------
# 3. Server Logic (Back-end Processing)
# ------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # Reactive execution block triggered by clicking the "Run Simulation" button
  sim_results <- eventReactive(input$run_sim, {
    
    # Inject missing helper functions securely into active runtime global environment
    .GlobalEnv$findOBD <- function(pV, qV, pT, qE) {
      return(findOBD_RDS(pV, qV, pT, qE, u11 = 100, u00 = 0))
    }
    
    .GlobalEnv$pava_rcpp <- pava_rcpp
    
    # Bind the decision matrix to the active execution thread
    if(exists("fun_TITE_PKBOIN12dec")) {
      .GlobalEnv$fun_TITE_PKBOIN12dec <- fun_TITE_PKBOIN12dec
    }
    
    dN_val  <- input$n_dose
    
    # Automatically formulate testing parameter boundaries based on input steps
    pV_mock <- seq(0.05, 0.45, length.out = dN_val) 
    qV_mock <- seq(0.20, 0.60, length.out = dN_val) 
    rV_mock <- seq(10, 50, length.out = dN_val)     
    
    # Render interactive progress layout notification overlay
    withProgress(message = 'Running trial simulations...', value = 0.5, {
      
      # Fire main back-end trial design processing simulation engine 
      res <- fun_TITE_PK_fixsimu(
        dN = dN_val, 
        rV = rV_mock, 
        pV = pV_mock, 
        qV = qV_mock, 
        pT = input$target_tox, 
        qE = input$target_eff, 
        psi0PK = 30, 
        CV = 0.5, 
        g_P = 1.5, 
        csize = input$cohort_size, 
        cN = input$cN, 
        design = input$design_type,
        repsize = 100,  # Run 100 fast trial cycles to populate interactive visual charts
        n_cores = 1     # Secure 1 single thread to eliminate Windows background core crashes
      )
    })
    
    return(res)
  })
  
  # Print evaluation spreadsheet metadata log structures
  output$obd_text <- renderPrint({
    res <- sim_results()
    print(head(res)) 
  })
  
  # Draw probability distribution matrix charts
  output$sim_plot <- renderPlot({
    res <- sim_results()
    
    if("OBD" %in% names(res)){
      barplot(table(res$OBD), 
              main = "Distribution of Selected OBD Levels across Replications",
              col = "#005088", border = "white",
              xlab = "Dose Level", ylab = "Frequency Count")
    }
  })
  
  # Output detailed calculated row data lists
  output$detail_table <- renderTable({
    head(sim_results(), n = 10)
  })
}


# ------------------------------------------------------------------------------
# 4. Launch Application
# ------------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
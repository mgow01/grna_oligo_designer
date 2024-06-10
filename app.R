library(shiny)
library(dplyr)
library(stringr)
library(Biostrings)

ui <- fluidPage(
  titlePanel("gRNA oligo designer for gRNAs to be ligated into BsaI-digested Cas9-YFP vector (pG474)"),
  p("Please paste the gR number(s) and gRNA sequence(s) in the field below"),
  p("Each gRNA should be on a separate row, with the gR number and sequence separated by a space or tab. As nucleic acids, only the characters ATCG are supported."),
  p("For example:"),
  p(em("gR1 TGAGACGCGCCGTCACGGCA")),
  p(em("gR2 GGCAGTGAGACGCGCCGTCA")),
  textAreaInput("entry", "",
                width = "100%",
                rows = 6),
  
  actionButton("submit", "Submit"),
  conditionalPanel(condition = "input.submit != 0",
                   textOutput("error_msg")),
  tableOutput("table"),
  uiOutput("downloadData"),
)

# Define server logic
server <- function(input, output, session) {
  rv <- reactiveValues(
    df = data.frame("gR Number" = as.numeric(),
                    "gRNA Sequence" = as.character(),
                    "Oligo Sequence to Order" = as.character(),
                    "Oligo Name" = as.character()
    ),
    seq_in = NULL
  )
  
  observeEvent(input$submit, {
    rv$seq_in <- unlist(str_extract_all(input$entry, "(?<=[[:blank:]])\\w+")) %>% 
      str_split(pattern = "") %>% unlist() %>% toupper()
    
    if (isTRUE(rv$seq_in %in% DNA_BASES %>% all())) {
      fw5 <- DNAString("AAGTT")
      fw3 <- DNAString("G")
      rv5 <- DNAString("AAAAC")
      rv3 <- DNAString("A")
      gR <- unlist(str_extract_all(input$entry, "gR\\d+"))
      gRNA <- unlist(str_extract_all(input$entry, "(?<=[[:blank:]])\\w+")) %>% 
        toupper()
      
      gRNA <- DNAStringSet(gRNA)
      gRNA_rc <- reverseComplement(gRNA)
      gRNA <- c(rbind(gRNA, gRNA_rc))
      gRNA <- DNAStringSet(gRNA)
      
      oligo <- gRNA
      for (i in seq_along(oligo)) {
        if (i %% 2 == 0) {
          oligo[[i]] <- c(rv5, oligo[[i]], rv3)
        }
        else {
          oligo[[i]] <- c(fw5, oligo[[i]], fw3)
        }
      }
      
      # labelling oligo names to indicate forward and reverse
      oligo_name <- rep(gR, each = 2)
      for (i in 1:length(oligo_name)) {
        if (i %% 2 == 0) {
          oligo_name[i] <- paste(oligo_name[i], "Rv")
        }
        else {
          oligo_name[i] <- paste(oligo_name[i], "Fw")
        }}
      
      #duplicating each gR number and then adding "revcom to the 2nd"
      gR <- rep(gR, each = 2)
      for (i in 1:length(gR)) {
        if (i %% 2 == 0) {
          gR[i] <- paste(gR[i], "reverse complemented")
        }
        else {
          gR[i] <- paste(gR[i])
        }}
      
      rv$df <- data.frame("gR Number" = gR,
                          "gRNA Sequence" = gRNA,
                          "Oligo Sequence to Order" = oligo,
                          "Oligo Name" = oligo_name,
                          check.names = FALSE)
    }
  })
  
  # output$error_msg <- (renderText({
  #   validate(need(isFALSE(rv$seq_in %in% DNA_BASES %>% all()),
  #                 message = "Error: input contains invalid DNA base"))
  # })
  # )
  
  output$error_msg <- renderText({
    if (isFALSE(rv$seq_in %in% DNA_BASES %>% all())) {
      print("Error: input contains invalid DNA base")
    } else {
      return()
    }
  })
  
  output$table <- renderTable({
    if (nrow(rv$df) < 1) {return()}
    else {rv$df}
  }
  )
  
  output$downloadData <- renderUI({
    if (nrow(rv$df) < 1) {return()}
    else {downloadButton("download", label = "Download table as csv")}
  })
  
  output$download <- downloadHandler(
    filename = function() {
      paste("gRNA oligo designs.csv")
    },
    content = function(file) {
      write.csv(rv$df, file, row.names = FALSE)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
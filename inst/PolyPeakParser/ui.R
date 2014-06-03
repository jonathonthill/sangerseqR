
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

shinyUI(pageWithSidebar(
  headerPanel(title = "", windowTitle = "Poly Peak Parser"),
  sidebarPanel(
    tags$head( tags$link(rel="stylesheet", type="text/css", href="styleupdates.css") ),
    tags$img(src="logo3.png", alt="Double Peak Parser", width="80%", class="unframed"),
    tags$hr(),
    tags$h4('1. Upload Data:'),
    conditionalPanel(condition = "!input.example", fileInput('seq', 'Upload Chromatogram File (.ab1 or .scf)')),
    checkboxInput('example', tags$span(style="color: red", 'Load Example Data'), FALSE),
    tags$h4('2. Set Chromatogram Options:'),
    helpText(tags$small("Note: Chromatogram will update live to reflect changes.")),
    sliderInput("x", strong("Approx. number of bases per row"), 10, 200, 100, 10, format=" ", ticks=FALSE),
    helpText(tags$small("Higher numbers fit more peaks on a single line.")),
    tags$div(class="float-left", numericInput("trim5", strong("5' Trim"), '30')),
    tags$div(class="float-left", numericInput("trim3", strong("3' Trim"), '30')),
    tags$br(style="clear:both"),
    checkboxInput('showtrim', 'Show Trimmed Region', FALSE),
    helpText(tags$small("Removes low quality bases from the ends.")),
    numericInput("ratio", strong('Signal Ratio Cutoff'), '0.33', min=0, max=1, step=.01),
    helpText(tags$small("Signal Ratio is calculated as peak signal/max signal for each position. Signals above this ratio are called as alternate bases.")),
    
    tags$h4('3. Enter Reference Sequence:'),
    conditionalPanel(condition = "!input.example", tags$textarea(id="ref", rows=8, cols=60,"")),
    conditionalPanel(condition = "input.example", tags$textarea(id="ref", rows=8, cols=60, "TGGCAAGAGAGCGACAGTCAGTCGGACTTACGAGTTGTTTTTACAGGCGCAATTCTTTTTTTAGAATATTATACATTCATCTGGCTTTTTGGATGCACCGATGAGAGATCCAGTTTTCACAGCGAACGCTATGGCTTATCACCCTTTTCACGCGCACAGGCCGGCCGACTTTCCCATGTCAGCTTTCCTTGCGGCGGCTCAACCTTCGTTCTTTCCAGCGCTCACTTTACCACCGGGTCTCAGTAAACCGCTGGCGGATCATGCGCTCTCCGGTGCGGCTGAAGCTGGTTTACACGCGGCGCTTGGACATCACCACCAGGCGGCTCATCTGCGCTCTTTCAAGGGTCTCGAGCCAGAGGAGGATGTTGAGGACGATCCTAAAGTTACATTAGAAGCTAAGGAGCTTTGGGATCAATTCCACAAAATTGGAACAGAAATGGTCATCACTAAATCAGGAAGGTAAGGTCTTTACATTATTTAACCTATTGAATGCTGCATAGGGTGATGTTATTATATTACTCTGCGAAGAGTTGGGTCTATTTTATCGTAAAATATACTTTACATTATAAAATATTGCTCGGTTAAAATTCAGATGTACTGGATGCTGACATAGCATCGAAGCCTCT"), tags$p(style="color: red;", "Example Data: Click Results Tab to see results.")),
    checkboxInput('trimref', 'Trim alt allele to match chromatogram', TRUE),
    helpText(tags$small("Sequence should include at least the region covered by the sequence results. Non-DNA characters (e.g. numbers) will automatically be removed.")),
    tags$hr(),
    tags$small('Created by the Yost Lab', tags$a(href="http://yost.genetics.utah.edu", "yost.genetics.utah.edu", target="_blank")),
    tags$br(),
    tags$small('Hosted by', tags$a(href="http://www.rstudio.com", target="_blank", "Rstudio")),
    tags$br(),
    tags$small('Using sangerseqR version ', textOutput('ssRversion'))
  ),
  
  mainPanel(
    tabsetPanel(id="maintabset",
                tabPanel("Instructions", includeHTML("instructions.html")),
                
                tabPanel("Chromatogram", 
                         conditionalPanel(
                           condition = "!output.fileUploaded", 
                           tags$p("A Chromatogram will be shown here when a valid sequencing file has been uploaded.")
                         ), 
                         conditionalPanel(
                           condition = "output.fileUploaded", 
                           plotOutput('chromatogram', height="auto")
                         )
                ),
                
                tabPanel("Results", 
                         conditionalPanel(
                           condition = "input.ref == '' && !input.example", 
                           tags$p("The Alternate Allele and Alignment will be shown here after a data file and reference sequence have been entered.")
                         ), 
                         conditionalPanel(
                           condition = "(output.fileUploaded && input.ref != '') || input.example", 
                           downloadButton('downloadData', 'Download PDF Report'),
                           tags$h4("Alignment"), 
                           verbatimTextOutput('alignment'),
                           tags$h4("Reference Sequence"), 
                           verbatimTextOutput('refseq'), 
                           tags$h4("Alternate Allele"), 
                           verbatimTextOutput('altseq'),
                           tags$h4("Alignment Header"), 
                           verbatimTextOutput('header')
                           
                         )
                )
    )
  )
))
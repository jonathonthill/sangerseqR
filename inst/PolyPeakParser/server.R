library(shiny)
library(knitr)
library(sangerseqR)
source('polypeakfunctions.R')

examplefile <- "example.ab1"
exampleref <- "TGGCAAGAGAGCGACAGTCAGTCGGACTTACGAGTTGTTTTTACAGGCGCAATTCTTTTTTTAGAATATTATACATTCATCTGGCTTTTTGGATGCACCGATGAGAGATCCAGTTTTCACAGCGAACGCTATGGCTTATCACCCTTTTCACGCGCACAGGCCGGCCGACTTTCCCATGTCAGCTTTCCTTGCGGCGGCTCAACCTTCGTTCTTTCCAGCGCTCACTTTACCACCGGGTCTCAGTAAACCGCTGGCGGATCATGCGCTCTCCGGTGCGGCTGAAGCTGGTTTACACGCGGCGCTTGGACATCACCACCAGGCGGCTCATCTGCGCTCTTTCAAGGGTCTCGAGCCAGAGGAGGATGTTGAGGACGATCCTAAAGTTACATTAGAAGCTAAGGAGCTTTGGGATCAATTCCACAAAATTGGAACAGAAATGGTCATCACTAAATCAGGAAGGTAAGGTCTTTACATTATTTAACCTATTGAATGCTGCATAGGGTGATGTTATTATATTACTCTGCGAAGAGTTGGGTCTATTTTATCGTAAAATATACTTTACATTATAAAATATTGCTCGGTTAAAATTCAGATGTACTGGATGCTGACATAGCATCGAAGCCTCT"

log <- "access.log"

shinyServer(function(input, output, session) {
  inputdata <- reactive({ 
    if(input$example) {
      return(makeBaseCalls(readsangerseq(examplefile), input$ratio))
    } else if(!is.null(input$seq)) {
      return(makeBaseCalls(readsangerseq(input$seq$datapath), input$ratio))
    } else return(NULL)
  })
  refseq <- reactive(
    if(input$example) cleanstring(exampleref)
    else cleanstring(input$ref)
  )
  h <- reactive({
    if (!is.null(input$seq) | input$example) {
      figheight(inputdata(), input$trim5, input$trim3, width=input$x, showtrim=input$showtrim)
    }
  }) 
 
  output$fileUploaded <- reactive({return(!is.null(input$seq) | input$example)})
  output$chromatogram <- renderPlot(
    chromatogram(inputdata(), showcalls="both", trim3=input$trim3, 
                 trim5=input$trim5, width=input$x, showtrim=input$showtrim, showhets=TRUE, cex.base=2), height=reactive(h())
  )
  output$h <- reactive(h())  
  outputdata <- reactive(alignchromatogram(inputdata(), 
                                           trim=input$trimref, 
                                           refseq=refseq(), 
                                           trim5=input$trim5, 
                                           trim3=input$trim3, 
                                           block.width=80
                                           )
                         )
  output$refseq <- renderText(
    if (nchar(input$ref) > 0 | input$example) {outputdata()$refseq})
  output$altseq <- renderText(
    if (nchar(input$ref) > 0 | input$example) {outputdata()$altseq})
  output$alignment <- renderText(
    if (nchar(input$ref) > 0 | input$example) {
      gsub(pattern="^.*#={39}(.+?)#-{39}.*$",
           replacement="\\1",
           x=outputdata()$alignment)
      
  })
  
  output$header <- renderText(
    if (nchar(input$ref) > 0 | input$example) {
      gsub(pattern="(^.+)#\\n#\\n#={39}.+$",
           replacement="\\1",
           x=outputdata()$alignment)
    }
  )

  
  observe({
    c(input$ratio, input$trim5, input$trim3)
    if(!is.null(input$seq) | input$example == TRUE)
      updateTabsetPanel(session, "maintabset", selected = "Chromatogram")
  })
  
  observe({
    inputref <- input$ref
    inputseq <- input$seq
    if(inputref != "" & !is.null(inputseq)) {
      updateTabsetPanel(session, "maintabset", selected = "Results")
    }
  })
  
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  output$downloadData <- downloadHandler(
    filename = function() paste(input$seq$name, "_report.pdf", sep=""),
    content = function(con) {
      pdfname = knit2pdf(input="peakparser_report.Rnw")
      file.copy(pdfname, con)
    }
  )
  
  #logging
  observe({
    if(!is.null(inputdata())) {
      isolate({
        alog <- file(log, "a")
        cat(format(Sys.time(), "%Y-%b-%d_%H-%M-%S"), file=alog)
        cat("\t", file=alog)
        if(input$example == TRUE) cat("Opened Example", file=alog)
        else cat(input$seq$name, file=alog)
        cat("\n", file=alog)
        close(alog)
      })
    }
  })
  output$ssRversion <- renderText(as.character(packageVersion("sangerseqR")))
})



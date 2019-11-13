#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(qqman)
library(dplyr)
df_refGene_chrom_txStart_name<-unique(read.csv("refGene.txt",sep="\t") %>% select(chrom,txStart,name2))
highlightedSnps<-c("rs145677414","rs145663603","rs145942463","rs146065087","rs146104244","rs146453317","rs146458139")
# library(ggplot2)
df_imp_hg19_Allpop_CLCLP_tdt<-read.csv("imp_hg19_Allpop_CLCLP_tdt.txt")
# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("imp_hg19_Allpop_CLCLP_tdt"),

   # Sidebar with a slider input for number of bins 
   #flowLayout(
   verticalLayout(
      
     selectInput("Chromosome", "Chromosome", c("ALL"="ALL","1"="1","2"="2","3"="3","4"="4","5"="5","6"="6","7"="7","8"="8","9"="9","10"="10","11"="11","12"="12","13"="13","14"="14","15"="15","16"="16","17"="17","18"="18","19"="19","20"="20","21"="21","22"="22","23"="23")
     ),
     #actionButton("dataChromosome", "Update Plot"),
     numericInput("Location", "Location", 10000000,1,243019906),
     #actionButton("dataLocation", "Update Plot"),
     #NonNumericselectInput("Flank_Size", "Flank_Size", c("100"="100", "1000"="1000", "10000"="10000", "100000"="100000", "1000000"="1000000"),100000),
     numericInput("FlankSize", "FlankSize", 1000000,100,121000000,1000000),
     #actionButton("dataFlankSize", "Update Plot"),
     #NO actionButton("go", "Update Plot"),
     # Show a plot of the generated distribution
     #TTT
     #textInput("rsNumber","SNP RS Number","rs12345"),
     #textInput("geneSymbol","Gene","ABC3"),
     textInput("rsNumber","SNP RS Number",""),
     textInput("geneSymbol","Gene",""),
     actionButton("go", "Update Plot"),
     actionButton("reset", "Reset Form"),
        plotOutput("qqmanPlot")
      
   
   )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  #dataChromosome <- eventReactive(input$makemanhattanChromosome, {input$Chromosome})
  #dataLocation <- eventReactive(input$makemanhattanLocation, {input$Location;})
  #dataFlankSize <- eventReactive(input$makemanhattanFlankSize, {input$FlankSize})
  #data <- eventReactive(input$go, {input$FlankSize})
  #NO rv<-reactiveValues(data=list(input$Chromosome,input$Location,input$FlankSize))
  #NO rv<-reactiveValues(data=list(8,1000000,10000))
  #TTT
  data <- eventReactive(input$go, {list(input$Chromosome,input$Location,input$FlankSize,input$rsNumber,input$geneSymbol)})
  output$qqmanPlot <- renderPlot({
    showNotification("Plot Being Created..",duration=5)
    ###TTT WORKS101317
    ### observeEvent(input$show, {
    ###  showNotification("Reset Form?",action = a(href = "javascript:location.reload();", "Reload page"))
    ### })

      # generate bins based on input$bins from ui.R
      #x    <- faithful[, 2] 
     
     #bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      #hist(x, breaks = bins, col = 'darkgray', border = 'white')
     #manhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==input$CHR),annotatePval = 0.005, annotateTop = TRUE,ylim = c(0, 10), col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F)
     #WORKS manhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==input$Chromosome),annotatePval = 0.000005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F)
     #WORKS manhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==input$Chromosome),annotatePval = 0.000005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(300000,900000))
     #WORKS manhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==input$Chromosome),annotatePval = 0.000005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(input$Location-500000,input$Location+500000))
     #WORKS manhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==input$Chromosome),annotatePval = 0.000005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(input$Location-input$FlankSize,input$Location+input$FlankSize))
     #WORKSmanhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==input$Chromosome),annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(input$Location-input$FlankSize,input$Location+input$FlankSize))
     #WORKSmanhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==input$Chromosome),annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(input$Location-input$FlankSize,input$Location+input$FlankSize),highlight=highlightedSnps)
     #NOobserveEvent(input$go,{manhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==dataChromosome()),annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(dataLocation()-dataFlankSize(),dataLocation()+dataFlankSize()),highlight=highlightedSnps)})
    #NOobserveEvent(input$go,{manhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==rv$data[[1]]),annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(rv$data[[2]]-rv$data[[3]],rv$data[[2]]+rv$data[[3]]),highlight=highlightedSnps)})
    #TTT_NO_101317_This shows up in browser:_Error:missing values are not allowed in subscripted assignments of data frames
    #manhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==input$Chromosome),annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(input$Location-input$FlankSize,input$Location+input$FlankSize),highlight=highlightedSnps)
    #if(data()[[1]]=="ALL"){
    if(data()[[1]]=="ALL" && data()[[4]]=="" && data()[[5]]==""){
        manhattan(df_imp_hg19_Allpop_CLCLP_tdt,annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,highlight=highlightedSnps)
    }else if(data()[[4]]!=""){
      #Example RS number: 3,rs9834484,99707246,2.736e-05
      #Modded RS number: rs151095213_ref_allele_AA
      # 8,rs148002472,100611199,8.208e-07
      # 17,rs111336083,8946894,3.086e-07
      # 19,chr19:33524932,33524932,4.2e-07
      if(nrow(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(SNP==data()[[4]]))>0){
        df_snpChosen<-df_imp_hg19_Allpop_CLCLP_tdt %>% filter(SNP==data()[[4]])
        #WORKSmanhattan(df_snpChosen,annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(df_snpChosen$BP-100,df_snpChosen$BP+100),highlight=highlightedSnps)
        #WORKSbutWRONGlocation: manhattan(df_snpChosen,annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(data()[[2]]-data()[[3]],data()[[2]]+data()[[3]]),highlight=highlightedSnps)
        manhattan(df_snpChosen,annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = T, genomewideline = T,xlim=c(df_snpChosen$BP-100,df_snpChosen$BP+100),highlight=highlightedSnps)
        showNotification(paste0("SNP ",data()[[4]]," found at hg19:",df_snpChosen$CHR,":",df_snpChosen$BP," pValue ",df_snpChosen$P),duration=15)
      }else{
        #showNotification("SNP Not Found",action = a(href = "javascript:location.reload();", "Reload page"))
        showNotification(paste0("SNP ",data()[[4]]," Not Found"),duration=15)
        manhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==data()[[1]]),annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(data()[[2]]-data()[[3]],data()[[2]]+data()[[3]]),highlight=highlightedSnps)
      }
    }else if(data()[[5]]!=""){
      if(nrow(df_refGene_chrom_txStart_name %>% filter(name2==data()[[5]]))>0){
        df_geneChosen<-df_refGene_chrom_txStart_name %>% filter(name2==data()[[5]])
        #NOmanhattan(df_geneChosen,annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(df_geneChosen$txStart-100,df_geneChosen$txStart+100),highlight=highlightedSnps)
        #WORKSmanhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==substr(df_geneChosen$chrom,4,6)),annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(df_geneChosen$txStart-100,df_geneChosen$txStart+100),highlight=highlightedSnps)
        #WORKSbutWRONGlocation: manhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==substr(df_geneChosen$chrom,4,6)),annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(data()[[2]]-data()[[3]],data()[[2]]+data()[[3]]),highlight=highlightedSnps)
        manhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==substr(df_geneChosen$chrom,4,6)),annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = T, genomewideline = T,xlim=c(df_geneChosen$txStart-100,df_geneChosen$txStart+100),highlight=highlightedSnps)
        showNotification(paste0("Gene ",data()[[5]]," txStart found at hg19:",df_geneChosen$chrom," ",df_geneChosen$txStart),duration=15)
      }else{
        #showNotification("Gene Not Found",action = a(href = "javascript:location.reload();", "Reload page"))
        showNotification("Gene ",data()[[5]]," Not Found",duration=15)
        manhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==data()[[1]]),annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(data()[[2]]-data()[[3]],data()[[2]]+data()[[3]]),highlight=highlightedSnps)
      }
    }else{
      #WORKS, WRONG, but this plot will override previous ifClause: manhattan(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR==data()[[1]]),annotatePval = 0.00005, annotateTop = FALSE,ylim = c(0, 10), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F,xlim=c(data()[[2]]-data()[[3]],data()[[2]]+data()[[3]]),highlight=highlightedSnps)
    }
    #TTT 101317NO-resets everything everytime observe({
    #   input$show
    #   updateSelectInput(session, "Chromosome", selected = 8)
    #   updateNumericInput(session, "Location", value = 10000000)
    #   updateNumericInput(session, "FlankSize", value = 1000000)
    #   updateTextInput(session, "rsNumber", value = "")
    #   updateTextInput(session, "geneSymbol", value = "")
    # })
    observeEvent(input$reset,{
        updateSelectInput(session, "Chromosome", selected = 8)
        updateNumericInput(session, "Location", value = 10000000)
        updateNumericInput(session, "FlankSize", value = 1000000)
        updateTextInput(session, "rsNumber", value = "")
        updateTextInput(session, "geneSymbol", value = "")
      })
   })
}

# Run the application 
#NOshinyApp(ui = ui, server = server, session = session)
shinyApp(ui = ui, server = server)


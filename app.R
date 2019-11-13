library(shiny)
library(qqman)
library(dplyr)
library(plotly)
library(manhattanly)
#source("http://bioconductor.org/biocLite.R")
#library(rtracklayer)
#library(Gviz)
# In order to create a UCSC genome browser like ideogram, we want to use Bioconductor
# package Gviz as described here:
#   
#   https://www.bioconductor.org/help/course-materials/2012/BiocEurope2012/GvizEuropeanBioc2012.pdf
#   
# We install it like so:
#   
#   source("http://bioconductor.org/biocLite.R")
# biocLite("Gviz")
# library(Gviz)
# biocLite("rtracklayer")
# library(rtracklayer)
# 
#   ..and using a fix found here:
#   
#   https://support.bioconductor.org/p/92643/
# 
#   ..and implemented for hg19 with these four lines of code:
#   
# cytoBandIdeo_hg19 <- read.table("cytoBandIdeo.txt", header=F, sep="\t")
# colnames(cytoBandIdeo_hg19) <-c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain')
# ideoTrack <- IdeogramTrack(genome="hg19",chromosome = "chr2", bands=cytoBandIdeo_hg19)
# plotTracks(ideoTrack, from=143910627, to=143938620, showBandId = TRUE, cex.bands = 0.5)
# 
# ..that file, cytoBandIdeo.txt, was obtained, and gunzipped, from here:
#   wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBandIdeo.txt.gz
# 
# ideoTrack <- IdeogramTrack(genome="hg19",chromosome = "chr2", bands=cytoBandIdeo_hg19)
# plotTracks(ideoTrack, from=143910627, to=143938620, showBandId = TRUE, cex.bands = 0.5)

#WORKScytoBandIdeo_hg19 <- read.table("cytoBandIdeo.txt", header=F, sep="\t")
#WORKScolnames(cytoBandIdeo_hg19) <-c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain')
#ideoTrack <- IdeogramTrack(genome="hg19",chromosome = "chr2", bands=cytoBandIdeo_hg19)
#plotTracks(ideoTrack, from=143910627, to=143938620, showBandId = TRUE, cex.bands = 0.5)


df_refGene_chrom_txStart_name <-
  unique(read.csv("refGene.txt", sep = "\t") %>% select(chrom, txStart, name2))

highlightedSnps <-
  c(
    "rs145677414",
    "rs145663603",
    "rs145942463",
    "rs146065087",
    "rs146104244",
    "rs146453317",
    "rs146458139"
  )
#df_imp_hg19_Allpop_CLCLP_tdt <-  read.csv("imp_hg19_Allpop_CLCLP_tdt/imp_hg19_Allpop_CLCLP_tdt.txt")
df_imp_hg19_Allpop_CLCLP_tdt <-
  read.csv("imp_hg19_Allpop_CLCLP_tdt.txt")

ui <- fluidPage(
  titlePanel("imp_hg19_Allpop_CLCLP_tdt"),
  verticalLayout(
    selectInput(
      "Chromosome",
      "Chromosome",
      c(
        "1" = "1",
        "2" = "2",
        "3" = "3",
        "4" = "4",
        "5" = "5",
        "6" = "6",
        "7" = "7",
        "8" = "8",
        "9" = "9",
        "10" = "10",
        "11" = "11",
        "12" = "12",
        "13" = "13",
        "14" = "14",
        "15" = "15",
        "16" = "16",
        "17" = "17",
        "18" = "18",
        "19" = "19",
        "20" = "20",
        "21" = "21",
        "22" = "22",
        "23" = "23"
        # "ALL" = "ALL"
        # "23" = "23",
        # "ALL" = "ALL"
      )
    ),
    numericInput("Location", "Location", 100000000, 1, 243019906),
    numericInput("FlankSize", "FlankSize", 100000000, 100, 121000000, 1000000),
    textInput("rsNumber", "SNP RS Number OR Gene Symbol (i.e. rs145677414 NTN1)", ""),
    #textInput("geneSymbol", "Gene", ""),
    actionButton("go", "Update Plot"),
    #actionButton("reset", "Reset Form"),
    #actionButton("zoomin", "Zoom In"),
    plotlyOutput("plotlyPlot"),
    plotOutput("qqmanPlot")
    #WORKSplotOutput("ideogramPlot")
  )
)
server <- function(input, output, session) {
  zoomLevel=0
  rv<-reactiveValues(data=list(
          1,
          100000000,
          100000000,
          "",
          ""
        ))
  observeEvent(input$go,{rv$data<-list(
          input$Chromosome,
          input$Location,
          input$FlankSize,
          input$rsNumber,
          input$geneSymbol
        )})
  #WORKS101717
  # data <-
  #   eventReactive(c(input$go),
  #     { list(
  #       input$Chromosome,
  #       input$Location,
  #       input$FlankSize,
  #       input$rsNumber,
  #       input$geneSymbol
  #     )
  #   })

  #WORKS101717_about to change all instances of __data()__ into this: __rv$data()
  #WORKSoutput$plotlyPlot <- renderPlotly(manhattanly(HapMap))
  #WORKSoutput$plotlyPlot <- renderPlotly(manhattanly(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR == 22)))
  #WORKSoutput$plotlyPlot <- renderPlotly(manhattanly(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR == 22),snp="SNP"))
  #WORKSoutput$plotlyPlot <- renderPlotly(manhattanly(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR == rv$data[[1]]),snp="SNP"))
  output$qqmanPlot <- renderPlot({
    #showNotification(paste0("zoomLevel:",zoomLevel,"input$zoomin:",input$zoomin))
    #showNotification("Plot Being Created..", duration = 5)
    #if (rv$data[[1]] == "1" && rv$data[[4]] == "" &&
    #if (rv$data[[4]] == "" && rv$data[[5]] == "") {
    if (rv$data[[4]] == "") {
      showNotification(paste0("Plot Being Created hg19 ",rv$data[[1]],":",toString(rv$data[[2]] - rv$data[[3]]),"-",toString(rv$data[[2]] + rv$data[[3]])), duration = 5)
      output$plotlyPlot <- renderPlotly(manhattanly(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR == rv$data[[1]]),snp="SNP"))
      manhattan(
        df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR == rv$data[[1]]),
        annotatePval = 0.00005,
        annotateTop = FALSE,
        ylim = c(0, 10),
        cex = 0.6,
        cex.axis = 0.9,
        col = c("blue4", "orange3"),
        suggestiveline = F,
        genomewideline = F,
        xlim = c(rv$data[[2]] - rv$data[[3]], rv$data[[2]] + rv$data[[3]]),
        highlight = highlightedSnps
      )
      #WORKSideoTrack <- IdeogramTrack(genome="hg19",chromosome = rv$data[[1]], bands=cytoBandIdeo_hg19)
      #NOoutput$ideogramPlot <- plotTracks(ideoTrack, from=rv$data[[2]] - rv$data[[3]], to=rv$data[[2]] + rv$data[[3]], showBandId = TRUE, cex.bands = 0.5)
      #NOoutput$ideogramPlot <- plotTracks(ideoTrack, from=rv$data[[2]] - rv$data[[3]], to=rv$data[[2]] + rv$data[[3]], showBandId = TRUE, cex.bands = 0.5)[[1]]
      #NOoutput$ideogramPlot <- plotTracks(ideoTrack, from=rv$data[[2]] - rv$data[[3]], to=rv$data[[2]] + rv$data[[3]], showBandId = TRUE, cex.bands = 0.5, panel.only = T)
      
    } else if (rv$data[[4]] != "") {
      geneOrSnpChosen<-gsub("^\\s+|\\s+$", "", rv$data[[4]])
      if (stringr::str_detect(tolower(geneOrSnpChosen), "^rs")) {
        if (nrow(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(SNP == tolower(geneOrSnpChosen))) >
            0) {
          df_snpChosen <-
            df_imp_hg19_Allpop_CLCLP_tdt %>% filter(SNP == tolower(geneOrSnpChosen))
          showNotification("Plot Being Created..", duration = 5)
          output$plotlyPlot <- renderPlotly(manhattanly(df_snpChosen,snp="SNP"))
          manhattan(
            df_snpChosen,
            annotatePval = 0.00005,
            annotateTop = FALSE,
            ylim = c(0, 10),
            cex = 0.6,
            cex.axis = 0.9,
            col = c("blue4", "orange3"),
            suggestiveline = T,
            genomewideline = T,
            xlim = c(df_snpChosen$BP - 100, df_snpChosen$BP + 100),
            highlight = highlightedSnps
          )
          showNotification(
            paste0(
              "SNP ",
              rv$data[[4]],
              " found at hg19:",
              df_snpChosen$CHR,
              ":",
              df_snpChosen$BP,
              " pValue ",
              df_snpChosen$P
            ),
            duration = 15
          )
        } else{
          showNotification(paste0("SNP ", tolower(geneOrSnpChosen), " Not Found"), duration = 15)
          
          # NOTWORKING-here we tried to render a manhattanly but that won't work
          #  manhattanly is plotly and ploty can be used in a shiny app but only 
          #  with the renderPlotly (NOT renderPlot), and it fills up a UI plot widget
          #  named plotlyOutput (NOT plotOutput), and so we try again below the 
          #  renderPlot portion of the code to build the manhattanly in a 
          # if(geneOrSnpChosen=="rs777"){
          #   showNotification("We show manhattanly HapMap", duration = 15)
          #   #NO102317 manhattanly(HapMap,snp = "SNP", gene = "GENE",annotation1 = "ZSCORE", annotation2 = "EFFECTSIZE",highlight = significantSNP)
          #   manhattanly(HapMap,snp = "SNP", gene = "GENE",annotation1 = "ZSCORE", annotation2 = "EFFECTSIZE",highlight = significantSNP)
          #   #manhattan(df_imp_hg19_Allpop_CLCLP_tdt)
          # }
          
          showNotification(paste0("Plot Being Created hg19 ",rv$data[[1]],":",toString(rv$data[[2]] - rv$data[[3]]),"-",toString(rv$data[[2]] + rv$data[[3]])), duration = 5)
          output$plotlyPlot <- renderPlotly(manhattanly(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR == rv$data[[1]]),snp="SNP"))
          manhattan(
            df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR == rv$data[[1]]),
            annotatePval = 0.00005,
            annotateTop = FALSE,
            ylim = c(0, 10),
            cex = 0.6,
            cex.axis = 0.9,
            col = c("blue4", "orange3"),
            suggestiveline = F,
            genomewideline = F,
            xlim = c(rv$data[[2]] - rv$data[[3]], rv$data[[2]] + rv$data[[3]]),
            highlight = highlightedSnps
          )
        }
        
        
        # This is dicey because:
        #  toupper(geneOrSnpChosen) will match rs1234, 
        #   this is okay because we already checked if the first two letters were RS/rs
        #  we cannot match any gene that starts with RS
      } else if (stringr::str_detect(toupper(geneOrSnpChosen), "^[A-Z]")) {
        if (nrow(df_refGene_chrom_txStart_name %>% filter(name2 == toupper(geneOrSnpChosen))) >
            0) {
          df_geneChosen <-
            df_refGene_chrom_txStart_name %>% filter(name2 == toupper(geneOrSnpChosen))
          showNotification("Plot Being Created..", duration = 5)
          output$plotlyPlot <- renderPlotly(manhattanly(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR == substr(df_geneChosen$chrom, 4, 6)),snp="SNP"))
          manhattan(
            df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR == substr(df_geneChosen$chrom, 4, 6)),
            annotatePval = 0.00005,
            annotateTop = FALSE,
            ylim = c(0, 10),
            cex = 0.6,
            cex.axis = 0.9,
            col = c("blue4", "orange3"),
            suggestiveline = T,
            genomewideline = T,
            xlim = c(
              df_geneChosen$txStart - 100,
              df_geneChosen$txStart + 100
            ),
            highlight = highlightedSnps
          )
          showNotification(
            paste0(
              "Gene ",
              toupper(geneOrSnpChosen),
              " txStart found at hg19:",
              df_geneChosen$chrom,
              " ",
              df_geneChosen$txStart
            ),
            duration = 15
          )
        } else{
          showNotification(paste0("Gene ", toupper(geneOrSnpChosen), " Not Found"), duration = 15)
          showNotification(paste0("Plot Being Created hg19 ",rv$data[[1]],":",toString(rv$data[[2]] - rv$data[[3]]),"-",toString(rv$data[[2]] + rv$data[[3]])), duration = 5)
          output$plotlyPlot <- renderPlotly(manhattanly(df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR == rv$data[[1]]),snp="SNP"))
          manhattan(
            df_imp_hg19_Allpop_CLCLP_tdt %>% filter(CHR == rv$data[[1]]),
            annotatePval = 0.00005,
            annotateTop = FALSE,
            ylim = c(0, 10),
            cex = 0.6,
            cex.axis = 0.9,
            col = c("blue4", "orange3"),
            suggestiveline = F,
            genomewideline = F,
            xlim = c(rv$data[[2]] - rv$data[[3]], rv$data[[2]] + rv$data[[3]]),
            highlight = highlightedSnps
          )
        }
        updateTextInput(session, "rsNumber", value = "")
        #updateTextInput(session, "geneSymbol", value = "")
      }
      
      
    } else{
      showNotification("Should never get here, bug.", duration = 5)
    }
    observeEvent(input$reset, {
      updateSelectInput(session, "Chromosome", selected = 8)
      updateNumericInput(session, "Location", value = 10000000)
      updateNumericInput(session, "FlankSize", value = 1000000)
      updateTextInput(session, "rsNumber", value = "")
      updateTextInput(session, "geneSymbol", value = "")
    })
  })
}
shinyApp(ui = ui, server = server)

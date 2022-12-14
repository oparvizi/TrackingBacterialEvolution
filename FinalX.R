library(rio)
library(ape)
library(ggtree)
library(ggplot2)
library(lubridate)
library(tidyr)
library(dplyr)
library(ggmap)
library(RColorBrewer)
library(dygraphs)
library(xts)
library(leaflet)
library(ggtree)
library(shiny)
library(shinythemes)
library(shinydashboard)

##### CALCULATING
##--------------------------------------------------------------------
# PHYLOGENETIC TREE (NJ)

##library(DECIPHER)
##library(msa)
##dna<- readDNAStringSet("Seqwithout.fasta") #read fasta data
##AT <- AlignTranslation(dna, type="AAStringSet") # align the translation
#BrowseSeqs(AT, highlight=1) # view the alignment
#DNA <- AlignSeqs(dna) # align the sequences directly without translation
#writeXStringSet(DNA, file="<<path to output file>>") # write the aligned sequences to a FASTA file
##dnaAln <- msa(dna)
##library(seqinr)
##dnaAln2 <- msaConvert(dnaAln, type="seqinr::alignment") #convert to seqinr
##d <- dist.alignment(dnaAln2, "identity") # calculate the distance 
##library(ape)
##myTree <- nj(d)
##tree<-ggtree(myTree, right = TRUE);tree
##saveRDS(file="./data/root.Seqwithout.RDS",ggtree(myTree))
##saveRDS(file="./data/circ.Seqwithout.tree.RDS",ggtree(myTree,layout="circular"))
##saveRDS(file="./data/dl.Seqwithout.tree.RDS",ggtree(myTree,layout="daylight"))
##saveRDS(file="./data/inCirc.Seqwithout.tree.RDS",ggtree(myTree,llayout="inward_circular"))

# Load data------------------------------------------------------
library(rio)
metadata_campy <- import("Sweden_Campy_metadata.xlsx")
export(metadata_campy, "Sweden_Campy_metadata.RDS")
metadata_campy<-readRDS("Sweden_Campy_metadata.RDS")

#Define UI-------------------------------------------------------
require(shiny)
require(shinydashboard)
require(leaflet)
require(dygraphs)
ui<-dashboardPage(skin = "red",
                  dashboardHeader(title = "Shiny Dashboard", titleWidth =200),
                  dashboardSidebar(
                    width = 200,
                    h3("Tree Options"),
                    radioButtons(inputId="treeLayout","Layout",
                                 choices=c("Rectanular"="rec",
                                           "Circular"="circ"),
                                 selected="rec"),
                    selectizeInput(inputId="colorBy",
                                   label="Color By",
                                   choices=c("region", "source"),
                                   multiple=FALSE,
                                   selected="region")
                  ),
                  dashboardBody(
                    tags$head(
                      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
                    ),
                    box(title="Phylogenic Tree",
                        width=6,
                        plotOutput("treePlot",
                                   brush = "plot_brush")#,
                        #verbatimTextOutput("info")
                    ),
                    box(title="Geographic Coordinate",
                        leafletOutput("caseMap"),
                        width=6),
                    
                    box(title="Timeline",
                        dygraphOutput("timeline"),
                        width=12),
                  )
)

#Define Server---------------------------------------------------
# This is the server logic for a Shiny web application.
library(ape)
library(ggtree)
library(lubridate)
library(tidyr)
library(dplyr)
library(ggmap)
library(RColorBrewer)
library(dygraphs)
library(xts)
library(leaflet)

server<-function(input, output) {
  
  #Function and variables------------------------------------------------------
  library(RColorBrewer)
  coul <- brewer.pal(4, "PuOr");coul
  ReSo<-metadata_campy %>% select(c("id","region","source" ));ReSo
  colR <- colorRampPalette(coul)(length(unique (metadata_campy$region)));colR
  colS <- colorRampPalette(coul)(length(unique (metadata_campy$source)));colS
  regionCol <- data.frame(colR,ord = unique (metadata_campy$region));regionCol
  sourceCol <- data.frame(colS,ord = unique (metadata_campy$source));sourceCol
  ReSo$region<-factor(ReSo$region, levels = regionCol$ord);ReSo$region
  ReSo$source<-factor(ReSo$source, levels = sourceCol$ord);ReSo$source
   
  colorTreeTip = function(tree,metadata_campy,var) {
     if(var %in% c("region")){
       t<-tree %<+% ReSo + geom_tippoint(mapping = aes(color = region), size=5, alpha=0.35) + 
             geom_tiplab(linetype = "dashed",align = TRUE) + #geom_text2(aes(label=tree$data$label)) +
             theme(legend.justification = c("right", "bottom")) +#, legend.text = element_text(size = 8))
             scale_color_manual(values=regionCol$colR, drop=FALSE,na.translate = FALSE);t
    }
    else if(var %in% c("source")){
      t<-tree %<+% ReSo + geom_tippoint(mapping = aes(color = source), size=5, alpha=0.35) + 
            geom_tiplab(linetype = "dashed",align = TRUE) + 
            theme(legend.justification = c("right", "bottom")) +#, legend.text = element_text(size = 8))
            scale_color_manual(values=sourceCol$colS, drop=FALSE,na.translate = FALSE);t
    }
    t + ggexpand(0.3, side = "h")
  }  
  ##### REACTIVE VARIABLES
  # metadata variable that changes reactive according to the
  # timeline date range
  metadataReactive <- reactive({
    startDate<-input$timeline_date_window[[1]]
    endDate<-input$timeline_date_window[[2]]
    
    if(is.null(startDate)){
      metadata_campy 
    }else{
      metadata_campy %>% filter(date>=startDate & date <= endDate)
    }
  })
  ##### VISUALIZATIONS
  ##--------------------------------------------------------------------
  # PHYLOGENETIC TREEE
  
  output$treePlot <- renderPlot({
    # Load trees that have already been stored.
    tree<-readRDS("./data/root.Seqwithout.RDS")  # default is rooted tree
    if(input$treeLayout=="circ"){
      tree<-readRDS("./data/circ.Seqwithout.tree.RDS") #alternative
    }
    
    tree<-colorTreeTip(tree,metadata_campy,input$colorBy)
    
    #Works - but buggy brushing interaction
    if(!is.null(input$plot_brush) & input$treeLayout == "rec"){
      e<-input$plot_brush
      tree<- tree + 
        xlim(e$xmin,e$xmax) +
        ylim(e$ymin,e$ymax)
    }
    #return the tree
    tree
  })
  
  # Little bit of testing code that shows what is being clicked on in the phylogenetic tree
  #
  output$info <- renderText({
    xy_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
    }
    xy_range_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1), 
             " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
    }
    
    paste0(
      "click: ", xy_str(input$plot_click),
      "dblclick: ", xy_str(input$plot_dblclick),
      "hover: ", xy_str(input$plot_hover),
      "brush: ", xy_range_str(input$plot_brush)
    )
  })
  
  ##### VISUALIZATIONS
  ##--------------------------------------------------------------------
  # MAP THAT SHOWS CASE COUNT
  
  output$caseMap<-renderLeaflet({
    m<-NULL
    metadata_campy<-readRDS("Sweden_Campy_metadata.RDS")
    if(input$colorBy=="region"){
      pal<-colorFactor(colorRampPalette(coul)(length(unique (metadata_campy$region))), domain = unique (metadata_campy$region)) #leaflet
      
      aggDat<-metadataReactive() %>%
        filter(region !="?") %>%
        group_by(region,region_lon,region_lat) %>%
        dplyr::count()%>% 
        mutate(popup=sprintf("%s = %d cases",region,n))
      
      m<-leaflet(aggDat) 
      
      m %>%
        addTiles()%>% 
        addCircleMarkers(
          lng=~region_lon,
          lat= ~region_lat,
          radius=~sqrt(n)*2,
          color = ~pal(region),
          stroke = FALSE, fillOpacity = 0.7,
          label=~as.character(popup),
          labelOptions = labelOptions(noHide = T)# delete the lable made noHide = F
        )
    }else if(input$colorBy=="source"){
      pal<-colorFactor(colorRampPalette(coul)(length(unique (metadata_campy$source))), domain = unique (metadata_campy$source)) #leaflet
      
      aggDat<-metadataReactive() %>%
        filter(source !="?") %>%
        group_by(source,region,region_lon,region_lat) %>%
        dplyr::count()%>% 
        mutate(popup=sprintf("%s (%s) = %d cases",region,source,n))
      
      m<-leaflet(aggDat) 
      
      m %>%
        addTiles()%>% 
        addCircleMarkers(
          lng=~region_lon,
          lat= ~region_lat,
          radius=~sqrt(n)*2,
          color = ~pal(source),
          stroke = FALSE, fillOpacity = 0.7,
          label=~as.character(popup),
          labelOptions = labelOptions(noHide = T)# delete the lable made noHide = F
        )
    }
  })
  ##### VISUALIZATIONS
  ##--------------------------------------------------------------------
  # TIMELINE
  
  output$timeline<-renderDygraph({
    ######
    # To create the dygraph, first generate a xts series for *each* of the countries.
    #count cases by date, it's also aggregatge by *month* so we're going to 
    #create a new time variable
  if(input$colorBy=="region"){
    library(xts)
    YearMonth<-NULL
    timeseriesData<-metadata_campy %>%
      mutate(YearMonth=ymd(sapply(yearMonth,function(x){paste(x,"01",sep="-")}))) %>% 
      group_by(YearMonth)%>% 
      dplyr::count(country) %>%
      complete(YearMonth,country) %>% #make sure that all dates are represented
      mutate(n=replace(n,is.na(n),0)) #turn NAs from above command in zeros
    
    #create an xts object
    xtsObj<-c()
    for(i in unique(timeseriesData$country)){
      temp<-timeseriesData %>%
        filter(country == i)
       xtsObj<-cbind(xtsObj,xts(temp$n, temp$YearMonth))
    }
    #name out object, so that it plots the time series correctly
    colnames(xtsObj)<-unique(timeseriesData$country)
    #now make the the dygraph (yay!)
      dygraph(xtsObj,height=600) %>% 
      dyOptions(stackedGraph = TRUE,colors = "#008000", pointSize = 5) %>% 
      dyRangeSelector(height = 20, strokeColor = "") %>%
      dyAxis("y", label = "Nr. of Isolates")
   }else if (input$colorBy=="source"){
    library(xts)
    YearMonth<-NULL
    timeseriesData<-metadata_campy %>%
      mutate(YearMonth=ymd(sapply(yearMonth,function(x){paste(x,"01",sep="-")}))) %>% 
      group_by(YearMonth)%>% 
      dplyr::count(source) %>%
      complete(YearMonth,source) %>% #make sure that all dates are represented
      mutate(n=replace(n,is.na(n),0)) #turn NAs from above command in zeros
    #create an xts object
    xtsObj<-c()
    for(i in unique(timeseriesData$source)){
      temp<-timeseriesData %>%
         filter(source == i)
       xtsObj<-cbind(xtsObj,xts(temp$n, temp$YearMonth))
    }
    #name out object, so that it plots the time series correctly
    colnames(xtsObj)<-unique(timeseriesData$source)
    #now make the the dygraph (yay!)
    dygraph(xtsObj,height=600) %>% 
      dyOptions(stackedGraph = TRUE,colors = sourceCol$colS, pointSize = 3) %>% 
      dyRangeSelector(height = 20, strokeColor = "") %>%
      dyLegend(width = 400) %>%
      dyAxis("y", label = "Source of Isolates")
    }
  })
}

#Shiny App--------------------------------------
shinyApp(ui, server)


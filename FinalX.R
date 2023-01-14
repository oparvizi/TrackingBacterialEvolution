#Set directory
setwd("F:/campyR");getwd()
#Install necessary laboratory
#if(!require(rio)) install.packages("rio", repos = "http://cran.us.r-project.org")
#if(!require(readxl)) install.packages("readxl", repos = "http://cran.us.r-project.org")
#if(!require(ape)) install.packages("ape", repos = "http://cran.us.r-project.org")
#if(!require(tidyr)) install.packages("tidyr", repos = "http://cran.us.r-project.org")
#if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org")
#if(!require(ggmap)) install.packages("ggmap", repos = "http://cran.us.r-project.org")
#if(!require(RColorBrewer)) install.packages("RColorBrewer", repos = "http://cran.us.r-project.org")
#if(!require(dygraphs)) install.packages("dygraphs", repos = "http://cran.us.r-project.org")
#if(!require(xts)) install.packages("xts", repos = "http://cran.us.r-project.org")
#if(!require(leaflet)) install.packages("leaflet", repos = "http://cran.us.r-project.org")
#if(!require(lubridate)) install.packages("lubridate", repos = "http://cran.us.r-project.org")
#if(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
#if(!require(plotly)) install.packages("plotly", repos = "http://cran.us.r-project.org")
#if(!require(shiny)) install.packages("shiny", repos = "http://cran.us.r-project.org")
#if(!require(shinythemes)) install.packages("shinythemes", repos = "http://cran.us.r-project.org")
#if(!require(shinydashboard)) install.packages("shinydashboard", repos = "http://cran.us.r-project.org")
#if(!require(reactable)) install.packages("reactable", repos = "http://cran.us.r-project.org")
#if(!require(jsonlite)) install.packages("jsonlite", repos = "http://cran.us.r-project.org")
#if(!require(jsonlite)) install.packages("jsonlite", repos = "http://cran.us.r-project.org")
#if(!require(dygraphs)) install.packages("dygraphs", repos = "http://cran.us.r-project.org")
#if(!require(dendextend)) install.packages("dendextend", repos = "http://cran.us.r-project.org")
#if(!require(seqinr)) install.packages("seqinr", repos = "http://cran.us.r-project.org")
#Some additional empty commands if you need a special package
#if(!require()) install.packages("", repos = "http://cran.us.r-project.org")
#if(!require()) install.packages("", repos = "http://cran.us.r-project.org")
#if(!require()) install.packages("", repos = "http://cran.us.r-project.org")
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("ggtree") 
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager") 
#BiocManager::install("babette") 
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager") 
#BiocManager::install("Biostrings") 
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager") 
#BiocManager::install("BiocGenerics")
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager") 
#BiocManager::install("DECIPHER")
#if (!require("BiocManager", quietly = TRUE))install.packages("BiocManager") BiocManager::install("")
#if (!require("BiocManager", quietly = TRUE))install.packages("BiocManager") BiocManager::install("")
#install_formats() 

library(rio)
library(readxl)
library(ape)
library(Biostrings) 
library(BiocGenerics)
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
library(plotly)
library(reactable)
library(dendextend)
library(shiny) 
library(shinythemes)
library(shinydashboard)


##### Persistent data storage methods in Shiny apps
###--------------------------------------------------------------------
#   https://shiny.rstudio.com/articles/persistent-data-storage.html


##### Generate Tree
###--------------------------------------------------------------------
#
#obtain NJ and UPGMA
library(DECIPHER)
library(msa)
fasta_file<-"Seqwithout.fasta"
stopifnot(file.exists(fasta_file))
dna<- readDNAStringSet(fasta_file) #read fasta data
AT <- AlignTranslation(dna, type="AAStringSet") # align the translation
#BrowseSeqs(AT, highlight=1) # view the alignment
#DNA <- AlignSeqs(dna) # align the sequences directly without translation
#writeXStringSet(DNA, file="<<path to output file>>") # write the aligned sequences to a FASTA file
dnaAln <- msa(dna)
library(seqinr)
dnaAln2 <- msaConvert(dnaAln, type="seqinr::alignment") #convert to seqinr
d <- dist.alignment(dnaAln2, "identity") # calculate the distance 
library(ape)
my_nj <- ape::nj(d); #tree_nj<-ggtree(my_nj, right = TRUE);tree_nj
saveRDS(file="./data/root.nj.RDS",ggtree(my_nj))
saveRDS(file="./data/circ.nj.tree.RDS",ggtree(my_nj,layout="circular"))
#saveRDS(file="./data/sla.nj.tree.RDS",ggtree(my_nj,layout="slanted"))
#saveRDS(file="./data/eqA.nj.tree.RDS",ggtree(my_nj,layout="equal_angle"))
#saveRDS(file="./data/dl.nj.tree.RDS",ggtree(my_nj,layout="daylight"))

library(phangorn)
my_upgma <- phangorn::upgma(d); #tree_upgma<-ggtree(my_upgma, right = TRUE);tree_upgma
saveRDS(file="./data/root.upgma.RDS",ggtree(my_upgma))
saveRDS(file="./data/circ.upgma.tree.RDS",ggtree(my_upgma,layout="circular"))

#obtain a BEAST2 posterior with babette
#library(seqinr)
#fasta_file2<-"Seqwithout.fasta"
#stopifnot(file.exists(fasta_file2))
#fasta_file2 <- read.fasta(file=fasta_file2, as.string = TRUE )
#library(babette)
#tree_beast <- bbt_run(fasta_file2)
#plot_densitree(
#  output$alignment_trees,
#  alpha = 0.01,
#  consensus = as.character(c(1:4)),
#  cex = 2.0,
#  scaleX = TRUE,
#  scale.bar = FALSE
#)

# Load data------------------------------------------------------
#library(readxl)
library(rio)
#metadata_campy <- read_excel("Sweden_Campy_metadata.xlsx")
metadata_campy <- import("Sweden_Campy_metadata.xlsx")
saveRDS(metadata_campy, file="Sweden_Campy_metadata.RDS")
metadata_campy<-readRDS("Sweden_Campy_metadata.RDS")
metadata_campy

#Define UI-------------------------------------------------------
require(shiny)
require(shinydashboard)
require(leaflet)
require(dygraphs)
ui<- tagList(
  tags$style("html, body{background-color: black; overflow-x: hidden; overflow-y: auto;}
              .container{ width: 100%; heigh:auto; margin: 0 auto; padding: 0; }
              @media screen and (min-width: 800px){.container{ width: auto;}}"
  ),
  tags$div(class="container",
           
           # HEADER ------------------------------------------------------------------ 
           
           dashboardPage(
             skin = "red", 
             dashboardHeader(
               title = span(HTML(paste0("<p><span style='font-size:50;'>&#9763;</span>
                               <span style='font-size:35;'>HazardRadaR</span></p>")), "HazardRadaR"), 
               titleWidth =250),
             dashboardSidebar(width = 250,
                              sidebarMenu(
                                HTML(paste0(
                                  "<br>",
                                  "<img style = 'display: block; margin-left: auto;margin-right: auto;' src='https://www.sva.se/media/cyybfdr0/sva_logo_e.svg'; width = '150'; height='50'; background-color: white;>",
                                  "<p style = 'text-align: center;'><small> <a href='https://en.wikipedia.org/wiki/National_Veterinary_Institute_(Sweden)' 
          target='_blank'>National veterinary institute</a></small></p>",
                                  "<br>"
                                )),
                                menuItem(("Genomic Data Visualization"), tabName = "vision"), 
                                menuItem(("- Tree Types"), radioButtons(inputId="treeTypes","Types", choices=c("Neighbor-joining"="nj", "UPGMA"="upgma"),#, "Bayesian"="bayes"),
                                                                        selected="nj")),
                                menuItem(("- Tree Options"), radioButtons(inputId="treeLayout","Layout", choices=c("Rectanular"="rec", "Circular"="circ"), selected="rec")),
                                menuItem(("- Color Options"), radioButtons(inputId="colorBy","Color By", choices=c("Region"="region", "Source"="source"), selected="region")),
                                
                                menuItem("Data Query", tabName = "query", icon=icon("table") ),
                                menuItem("Statistics & Modeling", tabName = "statistics", icon = icon("stats", lib = "glyphicon")),
                                menuItem("Home", icon = icon("home"), href="https://www.sva.se/en/our-topics/research/research-projects-at-sva/foka/a-comprehensive-assessment-of-the-impact-of-campylobacter-positive-broilers-on-human-infection-from-farm-to-molecular-epidemiology/"),
                                menuItem( "FAQs", tabName = "faqs", icon = icon('question-circle')),# href="https://www.sva.se/1685?culture=en-US"),
                                menuItem("Releases", tabName = "releases", icon = icon("tasks")),
                                HTML(paste0(
                                  "<br>","<br>","<br>",
                                  "<table style='margin-left:auto; margin-right:auto;'>",
                                  "<tr>",
                                  "<td style='padding: 5px;'><a href='https://www.facebook.com/Statens.veterinarmedicinska.anstalt' target='_blank'><i class='fab fa-facebook-square fa-lg'>
         </i></a></td>",
                                  "<td style='padding: 5px;'><a href='https://www.linkedin.com/company/national-veterinary-institute-sweden' target='_blank'><i class='fab fa-linkedin fa-lg'>
         </i></a></td>",
                                  "<td style='padding: 5px;'><a href='https://twitter.com/SVAexpertmyndig' target='_blank'><i class='fab fa-twitter fa-lg'></i></a></td>",
                                  "<td style='padding: 5px;'><a href='https://www.https://www.instagram.com/sva/' target='_blank'><i class='fab fa-instagram fa-lg'></i></a></td>",
                                  "</tr>",
                                  "</table>"
                                ),
                                HTML(paste0(
                                  "<script>",
                                  "var today = new Date();",
                                  "var yyyy = today.getFullYear();",
                                  "</script>",
                                  "<p style = 'text-align: center;'><small>&copy; - <a href='https://github.com/oparvizi/TrackingBacterialEvolution' target='_blank'>HazardRadaR</a> - <script>document.write(yyyy);</script></small></p>")   
                                ))
                              )                     
             ),
             
             # BODY ------------------------------------------------------------------ 
             
             dashboardBody(
               tabItems(
                 tabItem(tabName = "vision",
                         fluidRow(
                           box(title="Timeline", dygraphOutput("timeline", width = "auto", height = "120"), width=12, height=200),
                         ),
                         fluidRow(
                           box(title="Phylogenic Tree",plotOutput("treePlot", width = "auto", height = "500"), width=6, height=600), #brush = "plot_brush"),
                           box(title="Geographic Coordinate", leafletOutput("caseMap", width = "auto", height = "500"), width=6, height=600),
                         ),
                         fluidRow(
                           box(title="Mirror Tree", plotOutput("mirror_tree",width = "auto", height = "500"), width=12, height=600),
                         ),
                         fluidRow(
                           box(title="Presence/Absence Genes", plotlyOutput("heatGenes",width = "auto", height = "150"), width=6, height=220),
                           box(title="Presence/Absence Traits", plotlyOutput("heatTraits",width = "auto", height = "150"), width=6, height=220),
                         )
                 ),
                 tabItem(tabName = "statistics",
                         fluidRow( 
                           box(title="Regions", width=6, height=500,
                               plotOutput("plot1")),# click = "plot_click")), 
                           box(title="Time Period", width=6, height=500,
                               plotOutput("plot2")),# click = "plot_click")), 
                         )
                 ),
                 tabItem(tabName = "query",
                         fluidRow(
                           box(title="Datasets", width=12, height=650,
                               reactableOutput("ddynamic")),#, style="overflow-x: scroll; overflow-y: scroll;"),
                         ),
                         fluidRow(
                           box(title="Genomic Datasets", width=12, height=450,
                               reactableOutput("gdynamic")),#, style="overflow-x: scroll; overflow-y: scroll;"),
                         )
                 ),
                 tabItem(tabName = "faqs",
                         includeMarkdown("www/help.md")
                 ),
                 tabItem(tabName = "releases",
                         includeMarkdown("www/releases.md")
                 )
               )
             )
           )
  )
)

#Define Server---------------------------------------------------

#Function and variables------------------------------------------------------
library(RColorBrewer)
coul <- brewer.pal(4, "PuOr");coul
ReSo<-metadata_campy %>% select(c("id","region","source" ));ReSo
regio <- factor(metadata_campy$region);regio
regiol <- levels(regio[1]);regiol
sour <- factor(metadata_campy$source);regio
sourl <- levels(sour[1]);sourl
colR <- colorRampPalette(coul)(length(unique (metadata_campy$region)));colR
colS <- colorRampPalette(coul)(length(unique (metadata_campy$source)));colS
regionCol <- data.frame(colR,ord = regiol);regionCol
sourceCol <- data.frame(colS,ord = sourl);sourceCol
ReSo$region<-factor(ReSo$region, levels = regionCol$ord);ReSo$region
ReSo$source<-factor(ReSo$source, levels = sourceCol$ord);ReSo$source
objR1=objR2=objR3=objR4=objR5=objR6=objR7=objR8=NULL
for (i in 1:length(metadata_campy$id)){
  if (metadata_campy$region[i] == regiol[1]){
    objR1[i] <- metadata_campy$id[i]
  }
  if (metadata_campy$region[i] == regiol[2]){
    objR2[i] <- metadata_campy$id[i]
  }
  if (metadata_campy$region[i] == regiol[3]){
    objR3[i] <- metadata_campy$id[i]
  }
  if (metadata_campy$region[i] == regiol[4]){
    objR4[i] <- metadata_campy$id[i]
  }
  if (metadata_campy$region[i] == regiol[5]){
    objR5[i] <- metadata_campy$id[i]
  }  
  if (metadata_campy$region[i] == regiol[6]){
    objR6[i] <- metadata_campy$id[i]
  }  
  if (metadata_campy$region[i] == regiol[7]){
    objR7[i] <- metadata_campy$id[i]
  }   
  if (metadata_campy$region[i] == regiol[8]){
    objR8[i] <- metadata_campy$id[i]
  } 
  
  x1 <- objR1[!is.na(objR1)]
  x2 <- objR2[!is.na(objR2)]
  x3 <- objR3[!is.na(objR3)]
  x4 <- objR4[!is.na(objR4)]
  x5 <- objR5[!is.na(objR5)]
  x6 <- objR6[!is.na(objR6)]
  x7 <- objR7[!is.na(objR7)]
  x8 <- objR8[!is.na(objR8)]
}

objS1=objS2=objS3=objS4=objS5=NULL

for (i in 1:length(metadata_campy$id)){
  if (metadata_campy$source[i] == sourl[1]){
    objS1[i] <- metadata_campy$id[i]
  }
  if (metadata_campy$source[i] == sourl[2]){
    objS2[i] <- metadata_campy$id[i]
  }
  if (metadata_campy$source[i] == sourl[3]){
    objS3[i] <- metadata_campy$id[i]
  }
  if (metadata_campy$source[i] == sourl[4]){
    objS4[i] <- metadata_campy$id[i]
  }
  if (metadata_campy$source[i] == sourl[5]){
    objS5[i] <- metadata_campy$id[i]
  }  
  
  
  y1 <- objS1[!is.na(objS1)]
  y2 <- objS2[!is.na(objS2)]
  y3 <- objS3[!is.na(objS3)]
  y4 <- objS4[!is.na(objS4)]
  y5 <- objS5[!is.na(objS5)]
}

colorTreeTip = function(tree,metadata_campy,var) {

  if(var %in% c("region")){

    
      t<-tree %<+% ReSo + 
        geom_tippoint(aes(color = region), size=5, alpha=1) + 
        geom_point2(aes(subset=(label %in% dput(as.character(x1)))), shape=21, size=5, fill='#E66101')+ # If you add more color, pay attention to the color level
        geom_point2(aes(subset=(label %in% dput(as.character(x2)))), shape=21, size=5, fill='#EF862B')+
        geom_point2(aes(subset=(label %in% dput(as.character(x3)))), shape=21, size=5, fill='#F9AB55')+
        geom_point2(aes(subset=(label %in% dput(as.character(x4)))), shape=21, size=5, fill='#E7B482')+
        geom_point2(aes(subset=(label %in% dput(as.character(x5)))), shape=21, size=5, fill='#C7AEB2')+
        geom_point2(aes(subset=(label %in% dput(as.character(x6)))), shape=21, size=5, fill='#A69BC9')+
        geom_point2(aes(subset=(label %in% dput(as.character(x7)))), shape=21, size=5, fill='#826BB1')+
        geom_point2(aes(subset=(label %in% dput(as.character(x8)))), shape=21, size=5, fill='#5E3C99')+
        geom_tiplab(linetype = "dashed",linesize = 0.5, align = TRUE) + 
        geom_treescale() +
        theme(legend.justification = c("right", "bottom")) +
        scale_color_manual(values=regionCol$colR, drop=FALSE)

  }
  else if(var %in% c("source")){
    
      t<-tree %<+% ReSo + 
        geom_tippoint(aes(color = source), size=5, alpha=1) +
        geom_point2(aes(subset=(label %in% dput(as.character(y1)))), shape=21, size=5, fill='#E66101')+ 
        geom_point2(aes(subset=(label %in% dput(as.character(y2)))), shape=21, size=5, fill='#F7A24A')+
        geom_point2(aes(subset=(label %in% dput(as.character(y3)))), shape=21, size=5, fill='#D7B19A')+
        geom_point2(aes(subset=(label %in% dput(as.character(y4)))), shape=21, size=5, fill='#9C8FC3')+
        geom_point2(aes(subset=(label %in% dput(as.character(y5)))), shape=21, size=5, fill='#5E3C99')+
        geom_tiplab(linetype = "dashed",linesize = 0.5, align = TRUE) + 
        geom_treescale() +
        theme(legend.justification = c("right", "bottom")) +
        scale_color_manual(values=sourceCol$colS, drop=FALSE)#;t
    }
  t + ggexpand(0.3, side = "h")
  
}

server<-function(input, output) {
  
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
  ##### CALCULATING
  ##--------------------------------------------------------------------
  # DATA QUERY
  
  output$ddynamic <- renderReactable({
    
    metadata_campy<-readRDS("Sweden_Campy_metadata.RDS")
    reactable(
      metadata_campy[1:length(metadata_campy$id), ],
      searchable = TRUE,
      filterable = TRUE,
      defaultPageSize = 5,
      paginationType = "simple",
      language = reactableLang(
        searchPlaceholder = "Search...",
        noData = "No entries found",
        pageInfo = "{rowStart} to {rowEnd} of {rows} entries",
        pagePrevious = "\u276e",
        pageNext = "\u276f",
        
        # Accessible labels for assistive technologies such as screen readers.
        # These are already set by default, but don't forget to update them when
        # changing visible text.
        pagePreviousLabel = "Previous page",
        pageNextLabel = "Next page"
      )
    )
  })
  output$gdynamic <- renderReactable({
    
    metadata_campy<-readRDS("Sweden_Campy_metadata.RDS")
    reactable(
      metadata_campy[1:length(metadata_campy$id), ],
      searchable = TRUE,
      filterable = TRUE,
      defaultPageSize = 2,
      paginationType = "simple",
      language = reactableLang(
        searchPlaceholder = "Search...",
        noData = "No entries found",
        pageInfo = "{rowStart} to {rowEnd} of {rows} entries",
        pagePrevious = "\u276e",
        pageNext = "\u276f",
        
        # Accessible labels for assistive technologies such as screen readers.
        # These are already set by default, but don't forget to update them when
        # changing visible text.
        pagePreviousLabel = "Previous page",
        pageNextLabel = "Next page"
      )
    )
  })
  ##### CALCULATING
  ##--------------------------------------------------------------------
  # STATISTICS
  
  output$plot1 <- renderPlot({
    library(dplyr)
    # Count number cases per country
    aggDat1<-metadata_campy %>%
      filter(country !="?") %>%
      group_by(region,region_lon,region_lat) %>%
      dplyr::count()%>%
      mutate(popup=sprintf("%s = %d cases",region,n)) #create a popup for the map
    aggDat1
    aggDat1$region <- as.character(aggDat1$region)
    # Here's a very quick look at what this command generates for us:
    library(hrbrthemes)
    library(kableExtra)
    options(knitr.table.format = "html")
    # Barplot
    aggDat1 %>%
      filter(!is.na(n)) %>%
      arrange(n) %>%
      tail(20) %>%
      mutate(region=factor(region, region)) %>%
      ggplot( aes(x=region, y=n) ) +
      geom_bar(stat="identity", fill="#69b3a2") +
      coord_flip() +
      theme_ipsum() +
      ggtitle("Ferequency of isolates in Region")
  })
  
  output$plot2 <- renderPlot({
    library(ggplot2)
    library(dplyr)
    library(hrbrthemes)
    # Count number cases per country 
    aggDat2<-metadata_campy %>%
      filter(country !="?") %>%
      group_by(date) %>%
      dplyr::count()%>%
      mutate(popup=sprintf("%s = %d cases",date,n)) #create a popup for the map
    aggDat2
    aggDat2$date <- as.Date(aggDat2$date)
    # Plot
    aggDat2 %>%
      tail(10) %>%
      ggplot( aes(x=date, y=n)) +
      geom_line( color="grey") +
      geom_point(shape=21, color="black", fill="#69b3a2", size=6) +
      theme_ipsum() +
      ggtitle("Evolution of Isolate in time period")
    #plot(aggDat2$year, aggDat2$n, ylab="Cases", xlab="Year", col="darkred", lwd= 5) 
    #lines(lowess(x=aggDat2$year, y=aggDat2$n, f=.5))
  }) 
  
  output$vision<-renderDataTable({})
  output$treeTypes<-renderDataTable({
    if(input$treeLayout=="nj"){
      
      ##### CALCULATING
      ##--------------------------------------------------------------------
      # PHYLOGENETIC TREE (NJ)
      
      library(DECIPHER)
      library(msa)
      dna<- readDNAStringSet("Seqwithout.fasta") #read fasta data
      AT <- AlignTranslation(dna, type="AAStringSet") # align the translation
      #BrowseSeqs(AT, highlight=1) # view the alignment
      #DNA <- AlignSeqs(dna) # align the sequences directly without translation
      #writeXStringSet(DNA, file="<<path to output file>>") # write the aligned sequences to a FASTA file
      dnaAln <- msa(dna)
      library(seqinr)
      dnaAln2 <- msaConvert(dnaAln, type="seqinr::alignment") #convert to seqinr
      d <- dist.alignment(dnaAln2, "identity") # calculate the distance 
      library(ape)
      my_nj <- ape::nj(d); #tree1<-ggtree(my_nj, right = TRUE);#tree
      saveRDS(file="./data/root.nj.RDS",ggtree(my_nj))
      saveRDS(file="./data/circ.nj.tree.RDS",ggtree(my_nj,layout="circular"))
      #saveRDS(file="./data/sla.nj.tree.RDS",ggtree(my_nj,layout="slanted"))
      #saveRDS(file="./data/eqA.nj.tree.RDS",ggtree(my_nj,layout="equal_angle"))
      #saveRDS(file="./data/dl.nj.tree.RDS",ggtree(my_nj,layout="daylight"))
      #saveRDS(file="./data/inCirc.nj.tree.RDS",ggtree(my_nj,layout="inward_circular"))
    }
    if(input$treeLayout=="upgma"){
      library(phangorn)
      my_upgma <- phangorn::upgma  ;#tree2<-ggtree(mymy_upgma_nj, right = TRUE);#tree
      saveRDS(file="./data/root.upgma.RDS",ggtree(my_upgma))
      saveRDS(file="./data/circ.upgma.tree.RDS",ggtree(my_upgma,layout="circular"))
      
    }else if(input$treeLayout=="bayes"){
    }
  })
  
  ##### VISUALIZATIONS
  ##--------------------------------------------------------------------
  # TIMELINE
  library(dygraphs)
  library(xts)
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
      #      dygraph(xtsObj,height=600) %>% 
      #        dyOptions(stackedGraph = TRUE,colors = "#008000", pointSize = 5) %>% 
      #        dyRangeSelector(height = 20, strokeColor = "") %>%
      #        dyAxis("y", label = "Nr. of Isolates")
      
      dygraph(xtsObj,height=600) %>%
        dyOptions(labelsUTC = TRUE, fillGraph=TRUE, fillAlpha=0.1, drawGrid = FALSE, colors="#D8AE5A") %>%
        dyRangeSelector(height = 20, strokeColor = "") %>%
        dyCrosshair(direction = "vertical") %>%
        dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.2, hideOnMouseOut = FALSE)  %>%
        dyRoller(rollPeriod = 1)
      
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
      dygraph(xtsObj) %>% 
        dyOptions(stackedGraph = TRUE,colors = sourceCol$colS, pointSize = 3) %>% 
        dyRangeSelector(height = 20, strokeColor = "") %>%
        dyLegend(width = 400) %>%
        dyAxis("y", label = "Source of Isolates")
    }
  })
  ##### VISUALIZATIONS
  ##--------------------------------------------------------------------
  # PHYLOGENETIC TREE
  
  library(ape)
  library(ggtree)
  library(lubridate)
  output$treePlot <- renderPlot({
    # Load trees that have already been stored.  
    if(input$treeTypes=="nj" & input$treeLayout=="rec" ){
      tree1<-readRDS("./data/root.nj.RDS")  # default is rooted tree
      tree<-colorTreeTip(tree1,metadata_campy,input$colorBy)
    }
    if(input$treeTypes=="nj" & input$treeLayout=="circ"){
      tree2<-readRDS("./data/circ.nj.tree.RDS") #alternative
      tree<-colorTreeTip(tree2,metadata_campy,input$colorBy)
    }
    if(input$treeTypes=="upgma" & input$treeLayout=="rec" ){
      tree3<-readRDS("./data/root.upgma.RDS")  #alternative
      tree<-colorTreeTip(tree3,metadata_campy,input$colorBy)
    }
    if(input$treeTypes=="upgma" & input$treeLayout=="circ"){
      tree4<-readRDS("./data/circ.upgma.tree.RDS") #alternative
      tree<-colorTreeTip(tree4,metadata_campy,input$colorBy)
    }
    
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
  
  library(leaflet)
  library(ggmap)
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
          stroke = FALSE, fillOpacity = 1,
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
          stroke = FALSE, fillOpacity = 1,
          label=~as.character(popup),
          labelOptions = labelOptions(noHide = T)# delete the lable made noHide = F
        )
    }
  })

  ##### VISUALIZATIONS
  ##--------------------------------------------------------------------
  # Mirror Tree  
  
  output$mirror_tree <- renderPlot({
    m=NULL
    fasta_file<-"Seqwithout.fasta"
    stopifnot(file.exists(fasta_file))
    dna<- readDNAStringSet(fasta_file) #read fasta data
    AT <- AlignTranslation(dna, type="AAStringSet") # align the translation
    #BrowseSeqs(AT, highlight=1) # view the alignment
    #DNA <- AlignSeqs(dna) # align the sequences directly without translation
    #writeXStringSet(DNA, file="<<path to output file>>") # write the aligned sequences to a FASTA file
    dnaAln <- msa(dna)
    library(seqinr)
    dnaAln2 <- msaConvert(dnaAln, type="seqinr::alignment") #convert to seqinr
    d <- dist.alignment(dnaAln2, "identity") # calculate the distance 
    
    library(dendextend)
    #the agglomeration method to be used. This should be (an unambiguous abbreviation of) 
    #one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), 
    #"mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
    
    # Make 2 dendrograms, using 2 different clustering methods
    m_t1 <- d %>% dist() %>% hclust( method="average" ) %>% as.dendrogram()
    m_t2 <- d %>% dist() %>% hclust( method="centroid" ) %>% as.dendrogram()
    # Custom these kendo, and place them in a list
    m<-max(length(m_t1),length(m_t2))
    if(m<=9){
      coul <- brewer.pal(9, "Set1")
      valu<-coul[1:m]
    }else{
      coul <- brewer.pal(9, "Set1")
      col<-colorRampPalette(coul)(m)
      valu<-coul[1:m]
    }
      dl <- dendlist(
        m_t1 %>% 
          set("labels_col", value = valu, k=m) %>%
          set("branches_lty", 1) %>%
          set("branches_k_color", value = valu, k = m),
        m_t2 %>% 
          set("labels_col", value = valu, k= m) %>%
          set("branches_lty", 1) %>%
          set("branches_k_color", value = valu, k = m)
      )
    # Plot them together
    tanglegram(dl, 
               common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
               margin_inner=7,
               lwd=2
    )
  })
  
  ##### VISUALIZATIONS
  ##--------------------------------------------------------------------
  # HESTMAP GENES/TRAITS
  
  output$heatGenes <- renderPlotly({
    library(shiny)
    library(heatmaply)
    library(plotly)
    library(dplyr)
    # Load data 
    metadata_campy<-readRDS("Sweden_Campy_metadata.RDS")
    meta<-metadata_campy %>% select(c("id","TetO", "GyrA", "srRNA_23","fluoroquinolone_genotypes_2",
                                      "macrolide_genotypes_1","macrolide_genotypes_2","tetracycline_genotypes_1", 
                                      "tetracycline_genotypes_2","tetracycline_phenotype"));
    meta[is.na(meta)] <- 0
    colnames(meta) <- gsub("\\.", " ", colnames(meta)) #;colnames(meta)
    
    # Matrix format
    mat1 <- meta
    rownames(mat1) <- (mat1 [,1])
    mat1  <- mat1  %>% dplyr::select(TetO, GyrA, srRNA_23)
    mat1  <- data.matrix(mat1 ); #mat 
    
    #heatmap(mat1, scale="column")
    heatmaply(mat1,
              dendrogram = "row",
              xlab = "Genes", ylab = "Isolates", 
              #main = "Variation of presence/absence genes",
              scale = "none",
              margins = c(50,0,0,10),
              grid_color = "white",
              grid_width = 0.01,
              titleX = FALSE,
              hide_colorbar = TRUE,
              branches_lwd = 0.2,
              label_names = c("Isolate:", "Genes/Trait:", "Value:"),
              fontsize_row = 7, fontsize_col = 7,
              labCol = colnames(mat1),
              labRow = rownames(mat1),
              heatmap_layers = theme(axis.line=element_blank()),
              cellnote = mat1,
              scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                low = "white", 
                mid = "gold", 
                high = "red", 
                midpoint = 1
              )
    )
    
  })
  output$heatTraits <- renderPlotly({
    library(shiny)
    library(heatmaply)
    library(plotly)
    library(dplyr)
    # Load data 
    metadata_campy<-readRDS("Sweden_Campy_metadata.RDS")
    meta<-metadata_campy %>% select(c("id","TetO", "GyrA", "srRNA_23","fluoroquinolone_genotypes_2",
                                      "macrolide_genotypes_1","macrolide_genotypes_2","tetracycline_genotypes_1", 
                                      "tetracycline_genotypes_2","tetracycline_phenotype"));
    meta[is.na(meta)] <- 0
    colnames(meta) <- gsub("\\.", " ", colnames(meta)) #;colnames(meta)
    
    # Matrix format
    mat2 <- meta
    rownames(mat2) <- (mat2 [,1])
    mat2  <- mat2  %>% dplyr::select(fluoroquinolone_genotypes_2,
                                     macrolide_genotypes_1, macrolide_genotypes_2, tetracycline_genotypes_1, 
                                     tetracycline_genotypes_2, tetracycline_phenotype)
    mat2  <- data.matrix(mat2); #mat 
    
    #heatmap(mat2, scale="column")
    heatmaply(mat2,
              dendrogram = "row",
              xlab = "Genes", ylab = "Traits", 
              #main = "Variation of presence/absence genes",
              scale = "none",
              margins = c(50,0,0,10),
              grid_color = "white",
              grid_width = 0.01,
              titleX = FALSE,
              hide_colorbar = TRUE,
              branches_lwd = 0.2,
              label_names = c("Isolate:", "Genes/Trait:", "Value:"),
              fontsize_row = 7, fontsize_col = 7,
              labCol = colnames(mat2),
              labRow = rownames(mat2),
              heatmap_layers = theme(axis.line=element_blank()),
              cellnote = mat2,
              scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                low = "white", 
                mid = "gold", 
                high = "red", 
                midpoint = 1
              )
    )
    
  })

}

#Shiny App--------------------------------------
shinyApp(ui, server)


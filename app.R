# Shiny app to encapsulate the workflow of computing overlap in homeranges of 
#   animals from a selected herd, and for a selected date range

library(shiny)
library(shinydashboard)
library(mapview)
library(sf)
library(sp)
library(RColorBrewer)
library(leaflet)
library(dplyr)
library(DT)
library(visNetwork)

#

source("helpers.R")

# need to handle the CRS for herds yet, current hardcode to OR
  or.prj <- 26911 # UTM11N NAD83
  or.crsprj <- CRS("+init=epsg:26911")  
  wa.crsprj <- CRS("+init=epsg:2285")  # WA state plane north is 2285, south is 2927
  wa.prj <- 2285
  wgs.crsproj <- CRS("+init=epsg:4326") # CRS used for GPS data files
  wgs.proj <- 4326
  id.prj <- 8826 # Idaho Transverse Mercator
  id.crsprj <- CRS("+init=epsg:8826") 

# pre-loaded gps data is here #
  load("data/appdata.rda")

# go ahead and add bio year to gps table #
  gps <- addBioYear(gps)
 
  # LPMS has a weirdo point so we screen out anything east of lon -113
  # Lostine has 12LO27 the famous wanderer that has no GPS location in Lostine but was marked as a lamb there, will cause errors
    gps <- gps %>% filter(Longitude < -113.0 & AnimalID != "12LO27")
  
  
# determine some bookends from the data #
  #   these will be needed to modify the button options
  herds <- unique(gps$Herd)
  avail.dates <- gps %>% group_by(Herd) %>% summarise(MinDate=min(acquisitiontime),MaxDate=max(acquisitiontime))
  avail.dates$MinDate <- strftime(avail.dates$MinDate, format="%Y-%m-%d", tz="UTC")
  avail.dates$MaxDate <- strftime(avail.dates$MaxDate, format="%Y-%m-%d", tz="UTC")
  
  herd.bioyear.list <- vector(mode="list", length=length(herds))
  names(herd.bioyear.list) <- herds
  for (i in 1:length(herds)) {
    list.name <- herds[i]
    herd.bioyear.list[[list.name]] <- subset(gps,Herd==list.name,select=BioYear) %>% unique()
  }
  
# Convert the dataframe to a spatial object. Note that the
# crs= 4326 parameter assigns a WGS84 coordinate system to the 
# spatial object
  gps.sf <- st_as_sf(gps, coords = c("Longitude", "Latitude"), crs = 4326)
  rm(gps) # free up some memory
  


  
ui <- dashboardPage(
  
      dashboardHeader(title = "Tri-State TnR Home Range Overlap Tool",titleWidth = 450),
      dashboardSidebar(
          img(src = "IMG_0368.JPG", height = 150, width = 200),
          #strong("Idaho, Oregon, Washington"),
          p("Compute home range or UD overlap between all animals in a given herd (BETA). Note: brownian bridge method is computationally intensive, have patience."),
          selectInput("selectHerd", label = "Bighorn Herd:", choices = herds,selected = ""),
          selectInput("selectKernel", label = "Kernel Function:", choices = c("Bivariate Normal","Brownian Bridge"),
                      selected = "Bivariate Normal"),
          sliderInput("contour.perc", label = "Contour percentage:", min = 0, 
                      max = 100, value = 50),
          
          selectInput("selectMetric", label = "Overlap metric:", choices = c("Area overlap","UD volume"),
                      selected = "Area overlap"),
          # selectInput("dataSelect", label = "Subset data:",
          #              choices = list("None"="all","Date range"="drange"), 
          #              selected = "all"),
         
          uiOutput("subsetSelect"),
                   
          actionButton("runAnalysis", "Compute Overlap")
      ),
  
      
  dashboardBody(
    
    tabsetPanel(type = "tabs",
                
                tabPanel("Home Ranges", leafletOutput("homemap",width="95%",height=800)),
                tabPanel("GPS Locations", leafletOutput("gpsmap",width="95%",height=800)),
                tabPanel("Overlap Matrix", plotlyOutput("plot", width="85%",height=800),
                         downloadButton("downloadPlot", "Download")),
                tabPanel("Overlap Network", textOutput("networklabel"),visNetworkOutput("networkplot", width="85%",height=800),
                         downloadButton("downloadNetworkPlot", "Download Plot")),
                tabPanel("Overlap Dendrogram", plotOutput("clusterdend", width="85%",height=800),
                         downloadButton("downloadDendPlot", "Download Plot"),
                         downloadButton("downloadClusterData", "Download Data")),
                tabPanel("Data Table", DT::dataTableOutput("datatable"),
                         downloadButton("downloadTableData", "Download")),
                tabPanel("Tool-Info",
                         #img(src = "IMG_0368.JPG", height = 300, width = 450),
                         br(),
                         hr(),
                         h4(strong("Tool Description")),
                         p(style="text-align: justify; font-size = 25px",
                           "This application will first compute home ranges for every GPS-collared animal in the selected herd using the 'adehabitatHR' package and either a bivariate normal or brownian bridge kernel function. Then, the amount of overlap between each animal (area or UD) is calculated and stored in a NxN matrix. The kernel functionand contour level used to 
    compute the home range from the utilization distribution are user-configurable using the inputs on the left. The first tab contains a map of the computed home ranges. The polygons are also attributed with some of the capture and health testing results by individual. 
    
    In the next tab, we plot the NxN matrix that contains the fraction of each animals home range (by row in the matrix) 
    contained in every other animal in the herd (by columns). This plot gets messy with large numbers of animals. But you can explore
    the data and look at patterns or the relationship between specific animals.
    
     "),
                         br(),
                         p(style="text-align: justify; font-size = 25px","Next, treating this NxN matrix as a weighted adjacency matrix we can use tools 
    from the 'igraph' network analysis package to model the network and identify clusters in the data by viewing this data as a 
    directed social network, with the connection between animals weighted by the amount of home range overlap. 
    Note that in a directed network,the connection A to B can be different than B to A, which matches our data. The adjacency 
    matrix community structure is computed using a hierarchical walktrap method and displayed as a 1) a network plot and 
    2) a dendrogram. Now that we have mapped out existing social groups based on overlap in individual 
    home ranges, we can add the results of testing for Movi to our display."),
                         br(),
                         p(style="text-align: justify; font-size = 25px","In the network plot in the 4th tab, PCR status 
    at capture is denoted by node shape, where 'detected' is a star, 
     'indeterminate' a triangle, 'not detected' a circle. ELISA status at capture is denoted by text color where red 
    is 'detected', yellow is 'indeterminate' and green is 'not detected'. The fill color corresponds to cluster membership,
    matching the dendrogram in the previous tab. This is an interactive plot.
    
    In the fourth tab with the dendrogram, note that the group colors in the network plot match 
    the leaf text color in the dendrogram."),
                         br(),
                         p(style="text-align: justify; font-size = 25px","About downloads: the tool allows for the download of the raw matrix in the last tab. The dendrogram tab allows for a download of the 
                           exact copy of the dendrogram plot along with the community membership data in .csv. The download of the matrix and plots are slightly different as  those shown
                           since we're implementing interactive plotting methods in those two tabs."),
                         tags$blockquote("This application is still under development. 
                            Email scott.peckham78@gmail.com for troubleshooting or bug reporting."),
                         hr()))
                
  )
)

server <- function(input, output) {
  
  # first deal with the reactive UI issue (i.e. we don't know number of years or dates)
  output$subsetSelect <- renderUI({
    
    #if (input$dataSelect == "drange"){
      mind=unlist(c(avail.dates[which(input$selectHerd %in% avail.dates$Herd),2]))
      maxd=unlist(c(avail.dates[which(input$selectHerd %in% avail.dates$Herd),3]))
      dateRangeInput("dates",label="Date Range for home range computation:",
                     start  = "2023-01-01",
                     end    = "2023-03-31",
                     min=mind,
                     max=maxd)
    #}
    
  })
  
  inputHerd <- reactive({
    input$selectHerd
  })
  # # determine appropriate projected coord system for function input#
  # getMapProj <- reactive({
  #   switch(input$selectHerd,
  #   "Burnt River"=or.prj,
  #   "Lookout Mountain" =or.prj,
  #   "Yakima Canyon" = wa.prj,
  #   "Cleman Mountain" = wa.prj)
  # 
  # })
  
  # # reactive GPS table
  # gpstableInput <- reactive({
  #   t1 <- as.POSIXct(input$dates[1],tz="UTC")
  #   t2 <- as.POSIXct(input$dates[2],tz="UTC")
  #   switch(input$dataSelect, 
  #          "all" = gps.sf,
  #          "drange" = subset(gps.sf, acquisitiontime>= t1 & acquisitiontime <= t2)) %>% filter(Herd %in% inputHerd())
  # })
  
  # execute when action button is clicked #
  # #########--------------------------------------
  # 
  observeEvent(input$runAnalysis, {
    cat("Computing home ranges for:", input$selectHerd, " for ",input$dates[1], " to ", input$dates[2])
  })
  
  getData <- eventReactive(input$runAnalysis, {
    t1 <- as.POSIXct(input$dates[1],tz="UTC")
    t2 <- as.POSIXct(input$dates[2],tz="UTC")
    outdata <- subset(gps.sf, acquisitiontime>= t1 & acquisitiontime <= t2) %>% filter(Herd %in% input$selectHerd)
    
  })
  
  analysisHR <- eventReactive(input$runAnalysis,{
    
    output.proj <-switch(input$selectHerd,
           "Burnt River"=or.prj,
           "Lookout Mountain" = or.prj,
           "Lostine" = or.prj,
           "Yakima Canyon" = wa.prj,
           "Cleman Mountain" = wa.prj,
           "Lower Salmon" = id.prj,
           "Lower Panther Main Salmon" = id.prj)
    crs.projection <- switch(input$selectHerd,
                             "Burnt River"= or.crsprj, "Lookout Mountain" = or.crsprj, "Lostine"= or.crsprj,
                             "Yakima Canyon" = wa.crsprj, "Cleman Mountain" = wa.crsprj,
                             "Lower Salmon" = id.crsprj,
                             "Lower Panther Main Salmon" = id.crsprj)
    
    withProgress(message = paste("Computing:", input$selectHerd, " over ",input$dates[1], " to ", input$dates[2]), 
                 value = .25, {
    # compute homeranges
      #homeranges <- suppressWarnings(calculateHomerange(getData(),min.fixes=30,contour.percent=input$contour.perc, output.proj=output.proj))
      if (input$selectKernel == "Bivariate Normal") {
        homeranges <- suppressWarnings(calculateHomerange(getData(),min.fixes=30,contour.percent=input$contour.perc, output.proj=output.proj))
       }
      if (input$selectKernel == "Brownian Bridge"){
        homeranges <- suppressWarnings(calculateBBHomerange(getData(),min.fixes=30,contour.percent=input$contour.perc, output.proj=output.proj))
      } 
      
      
      
    # Increment the progress bar, and update the detail text.
      incProgress(0.75, detail = "Computing overlap of home ranges...")
      if (input$selectMetric=="Area overlap"){
      overlap <- suppressWarnings(calculateHomerangeOverlap(homeranges$homeranges))
      }
      if (input$selectMetric=="UD volume"){
        overlap <- kerneloverlaphr(homeranges$ud, method = "PHR",
                        percent = input$contour.perc)
        overlap[row(overlap) == col(overlap)] <- NA   # the overlap function returns 1 on diagonal, NA is better for display
                                                      # doesn't affect calculation because I ignore diagonal in the igraph calls
      }
    
      
    })
    
    return(
      list(
        homeranges=homeranges$homeranges,
        overlap=overlap
      )
    )
  })
    
  overlapHR <- eventReactive(input$runAnalysis, {
    # compute overlap
    calculateHomerangeOverlap(analysisHR())
  })  
  
  overlapPlot <- function(){
    .tmp <-analysisHR()
    ov <- .tmp$overlap
    overlapImagePlot(ov)
  }
  
  MatrixPlotDownload <- function(){
    .tmp <-analysisHR()
    ov <- .tmp$overlap
    overlapMatrixPlotDownload(ov)
  }
  
  # ClusterPlot <- reactive({
  #   .tmp <-analysisHR()
  #   ov <- .tmp$overlap
  #   overlapClusterPlot(ov)
  # })
  ClusterDend <- function(){
    .tmp <-analysisHR()
    ov <- .tmp$overlap
    suppressWarnings(overlapClusterDend(ov))
  }
  
  NetworkPlot <- function(){
    .tmp <-analysisHR()
    ov <- .tmp$overlap
    suppressWarnings(overlapNetworkPlot(ov,gps.sf))
  }
 
  NetworkPlotDownload <- function(){
    .tmp <-analysisHR()
    ov <- .tmp$overlap
    suppressWarnings(overlapNetworkPlotOutput(ov,gps.sf))
  }
  
  ClusterData <- reactive({
    .tmp <-analysisHR()
    ov <- .tmp$overlap
    suppressWarnings(getClusterData(ov))
  })
  
  overlapTable <- reactive({
    .tmp <- analysisHR()
    .tmp$overlap
  })
  
  output$gpsmap <- renderLeaflet({
    mindate <- as.POSIXct(input$dates[1],tz="UTC")
    maxdate <- as.POSIXct(input$dates[2],tz="UTC")
    zcol <- "AnimalID"
    plotdata <- getData()
    n.animals <- length(unique(plotdata$AnimalID))
    n.herd <- length(unique(plotdata$Herd))

    #mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(n.animals)
    #
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    if(n.animals > length(col_vector)) mycolors <- sample(col_vector, n.animals, replace=TRUE) else mycolors <- sample(col_vector, n.animals)
    makeGPSMap(plotdata, zcol=zcol,colors=mycolors,alpha=0.8)
  })
  
  output$homemap <- renderLeaflet({
    .tmp <-analysisHR()
    homeR <- .tmp$homeranges
    
    #mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(homeR$id))
    # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    # mycolors <- sample(col_vector, length(homeR$id))
    # makeHomerangeMap(homeR,zcol="id",colors=mycolors,alpha=0.8)
    
    makeHomerangeMap(homeR)
    
  })
  
  output$plot <- renderPlotly({
    overlapPlot()
  })
  
  output$clusterdend <- renderPlot({
    ClusterDend()
  })

  output$networklabel <- renderText({
    "In the network plot below, PCR status at capture is denoted by node shape, where 'detected' is an octagon, 
     'indeterminate' a triangle, 'not detected' a circle. ELISA status at capture is denoted by text color where red 
    is 'detected', yellow is 'indeterminate' and green is 'not detected'. The fill color corresponds to cluster membership, matching the dendrogram in the previous tab."
    
    })
  
  # output$networkplot <- renderPlot({
  #   NetworkPlot()
  # })
  output$networkplot <- renderVisNetwork({
    NetworkPlot()
  })
  
  
  # render table using DT to allow scroll
  output$datatable <- DT::renderDataTable({
    .tmp <-analysisHR()
    DT::datatable(format(as.data.frame(.tmp$overlap),digits=2),options = list(
      pageLength=25, scrollX='400px'), filter = 'top')
  })
    
  # Downloadable matrix overlap plot ----
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste(input$selectHerd, "_OverlapMatrixPlot.pdf", sep = "")
    }, 
    content = function(file,width=13,height=10) {
      pdf(file)
      MatrixPlotDownload()
      dev.off()
    },
    contentType = "application/pdf"
  )
  
  # Downloadable network plot ----
  output$downloadNetworkPlot <- downloadHandler(
    filename = function() {
      paste(input$selectHerd, "_CommunityNetworkPlot.pdf", sep = "")
    },
    content = function(file) {
      pdf(file,width=12,height=10)
      NetworkPlotDownload()
      dev.off()
    },
    contentType = "application/pdf"
  )
  
  # Downloadable cluster plot ----
  output$downloadDendPlot <- downloadHandler(
    filename = function() {
      paste(input$selectHerd, "_ClusterDendrogram.pdf", sep = "")
    },
    content = function(file) {
      pdf(file,width=12,height=10)
      ClusterDend()
      dev.off()
    },
    contentType = "application/pdf"
  )

  # Downloadable csv of clusters ----
  output$downloadClusterData <- downloadHandler(
    filename = function() {
      paste(input$selectHerd,"HomerangeOverlapClusters.csv", sep = "")
    },
    content = function(file) {
      write.csv(ClusterData(), file, row.names = TRUE)
    }
  )
  
  # Downloadable csv of selected table dataset ----
  output$downloadTableData <- downloadHandler(
    filename = function() {
      paste(input$selectHerd,"HomerangeOverlapTable.csv", sep = "")
    },
    content = function(file) {
      write.csv(overlapTable(), file, row.names = TRUE)
    }
  )
 
}

shinyApp(ui, server,options = list(launch.browser = TRUE))

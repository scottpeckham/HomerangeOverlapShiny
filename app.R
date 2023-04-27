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

# pre-loaded gps data is here #
  load("data/appdata.rda")

# go ahead and add bio year to gps table #
  gps <- addBioYear(gps)
 
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
  gps.sf <- st_as_sf(gps, coords = c("longitude", "latitude"), crs = 4326)


  
ui <- dashboardPage(
  
      dashboardHeader(title = "Tri-State TnR Home Range Overlap Tool",titleWidth = 450),
      dashboardSidebar(
          img(src = "IMG_0368.JPG", height = 150, width = 200),
          strong("Idaho, Oregon, Washington"),
          p("Compute home range or UD overlap between all animals in a given herd (BETA). Note: brownian bridge method is computationally intensive, have patience."),
          #p("Note: brownian bridge method is computationally intensive, have patience."),
          selectInput("selectHerd", label = "Bighorn Herd:", choices = herds,selected = ""),
          selectInput("selectKernel", label = "Kernel Function:", choices = c("Bivariate Normal","Brownian Bridge"),
                      selected = "Bivariate Normal"),
          sliderInput("contour.perc", label = "Contour percentage:", min = 0, 
                      max = 100, value = 50),
          selectInput("selectMetric", label = "Overlap metric:", choices = c("Area overlap","UD volume"),
                      selected = "Area overlap"),
          selectInput("dataSelect", label = "Subset data:",
                       choices = list("None"="all","Date range"="drange"), 
                       selected = "all"),
         
          uiOutput("subsetSelect"),
                   
          actionButton("runAnalysis", "Compute Overlap")
      ),
  
      
  dashboardBody(
    
    tabsetPanel(type = "tabs",
                tabPanel("Home Ranges", mapviewOutput("homemap",width="95%",height=800)),
                tabPanel("GPS Locations", mapviewOutput("gpsmap",width="95%",height=800)),
                tabPanel("Overlap Matrix", plotlyOutput("plot", width="85%",height=800),
                         downloadButton("downloadPlot", "Download")),
                tabPanel("Overlap Dendrogram", plotOutput("clusterdend", width="85%",height=800),
                         downloadButton("downloadClusterPlot", "Download Plot"),
                         downloadButton("downloadClusterData", "Download Data")),
                tabPanel("Overlap Network", textOutput("networklabel"),visNetworkOutput("networkplot", width="85%",height=800),
                         downloadButton("downloadNetworkPlot", "Download Plot")),
                tabPanel("Data Table", DT::dataTableOutput("datatable"),
                         downloadButton("downloadTableData", "Download"))
                
              )
  )
)



server <- function(input, output) {
  
  # first deal with the reactive UI issue (i.e. we don't know number of years or dates)
  output$subsetSelect <- renderUI({
    
    if (input$dataSelect == "drange"){
      mind=unlist(c(avail.dates[which(input$selectHerd %in% avail.dates$Herd),2]))
      maxd=unlist(c(avail.dates[which(input$selectHerd %in% avail.dates$Herd),3]))
      dateRangeInput("dates",label="Date Range for home range computation:",
                     start  = "2023-02-10",
                     end    = "2023-03-20",
                     min=mind,
                     max=maxd)
    }
    
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
    switch(input$dataSelect, 
                     "all" = gps.sf,
                     "drange" = subset(gps.sf, acquisitiontime>= t1 & acquisitiontime <= t2)) %>% filter(Herd %in% input$selectHerd)
    
  })
  analysisHR <- eventReactive(input$runAnalysis,{
    
    output.proj <-switch(input$selectHerd,
           "Burnt River"=or.prj,
           "Lookout Mountain" =or.prj,
           "Yakima Canyon" = wa.prj,
           "Cleman Mountain" = wa.prj)
    crs.projection <- switch(input$selectHerd,
                             "Burnt River"=or.crsprj, "Lookout Mountain" =or.crsprj,
                             "Yakima Canyon" = wa.crsprj, "Cleman Mountain" = wa.crsprj)
    
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
      
      proj4string(homeranges$homeranges) <- crs.projection   # reset this, kernelBB screws up the projection for some reason
      
    # Increment the progress bar, and update the detail text.
      incProgress(0.75, detail = "Computing overlap of home ranges...")
      if (input$selectMetric=="Area overlap"){
      overlap <- suppressWarnings(calculateHomerangeOverlap(homeranges$homeranges))
      }
      if (input$selectMetric=="UD volume"){
        overlap <- kerneloverlaphr(homeranges$ud, method = "VI",
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
    mycolors <- sample(col_vector, n.animals)
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
    "In the plot of the network below, PCR capture status is denoted by shape where an octagon is detected, 
     triangle is indeterminate, and circle is not detected. ELISA status is denoted by text color where red 
    is detected, yellow is indeterminate and green is not detected. The fill color corresponds to cluster membership."
    
    })
  
  # output$networkplot <- renderPlot({
  #   NetworkPlot()
  # })
  output$networkplot <- renderVisNetwork({
    NetworkPlot()
  })
  # # typical table rendering
  # output$datatable <- renderTable({
  #   .tmp <-analysisHR()
  #   .tmp$overlap
  # },rownames=TRUE,striped=TRUE)

  
  # render table using DT to allow scroll
  output$datatable <- DT::renderDataTable({
    .tmp <-analysisHR()
    DT::datatable(format(as.data.frame(.tmp$overlap),digits=2),options = list(
      pageLength=25, scrollX='400px'), filter = 'top')
  })
    
  # Downloadable overlap plot ----
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste(input$selectHerd, "_HomerangeOverlapPlot.pdf", sep = "")
    }, 
    content = function(file) {
      pdf(file)
      overlapPlot()
      dev.off()
    },
    contentType = "application/pdf"
  )
  
  # Downloadable cluster plot ----
  output$downloadClusterPlot <- downloadHandler(
    filename = function() {
      paste(input$selectHerd, "_HomerangeClusterPlot.pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      ClusterPlot()
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

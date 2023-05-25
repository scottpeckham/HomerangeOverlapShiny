# Helper file for Homerange overlap Shiny app
# this script houses some functions to compute homeranges from GPS tables from our DB#
#   assumes lat/lon are in WGS84 (i.e. vectronics data format)

# required packages #
library(sf)
library(adehabitatHR)
library(adehabitatLT)
library(RColorBrewer)
library(mapview)
library(sp)
library(dplyr)
library(igraph)
library(heatmaply)
library(dendextend)
library(visNetwork)
library(colorRamps)
library(fields)
library(colorspace)
# required packages #

getBioYear <- function(indate){
  year <- as.numeric(strftime(indate, format="%Y",tz="UTC"))
  month <- as.numeric(strftime(indate, format = "%m",tz="UTC"))
  if(month > 4) bioyear <- year else bioyear <- year-1
  return(bioyear)
}
addBioYear <- function(gps.data){
  # create a capture year column #
  gps.data$BioYear <-unlist(lapply(gps.data$acquisitiontime, getBioYear))
  return(gps.data)
}

calculateHomerange <- function(gps, min.fixes, contour.percent=95, output.proj, output.UD=FALSE){
  
  # check number of locations per animal, adehabitat requires >5/per
    loc.per.animal <- gps %>% group_by(AnimalID) %>%
      summarize(Fix_count =  n())
    drops <- which(loc.per.animal$Fix_count < min.fixes,arr.ind=TRUE)
    excludes <- loc.per.animal$AnimalID[drops]
    if (length(drops)>0) gps <- gps[!gps$AnimalID %in% excludes,]
  
  # ensure AnimalID is a factor
    gps$AnimalID <- factor(gps$AnimalID)
  
  # create a table that holds attributes for each animal (to add back to polygons)
    att.table <- gps[match(unique(gps$AnimalID),gps$AnimalID),]
  
  gps.sf <- subset(gps,select= AnimalID) # Keep just the "AnimalID" column of spatial object (again a restraint of the package) 
 
   # Transform input gps to proj.crs (typically UTM for Oregon)
    gps.sf.proj <- st_transform(gps.sf,crs=output.proj)  # transform to user input or UTM11 WGS84
  
  # convert from sf to spatial points data frame (reqd by adehabitatHR)
    gps.collarIDT <- as(gps.sf.proj, "Spatial")
    
    print(paste("Kernel function: Bivariate normal"))
    print(paste("Total number of rows in GPS table for HR calculation: ",nrow(gps.sf)))
    print(paste("Date range for HR calculation: ",(range(gps$acquisitiontime))[1]," to ",(range(gps$acquisitiontime))[2]))
    print(paste("Contour level is set to: ",contour.percent, "%"))
  
  # transform each data set in list from WGS84 to projected coords
  # and compute UD and polygon HR
    kud <- kernelUD(gps.collarIDT, grid=400,same4all = TRUE, h="href",extent=0.5)
    homeranges <- getverticeshr(kud, percent=contour.percent)
  
  # add some attributes to the polygon data
    newtable <- left_join(homeranges@data,att.table,join_by(id == AnimalID))
    homeranges@data <- newtable[,1:10]
    
  output <- list(homeranges=homeranges,ud=kud)
  return(output)
  
}

calculateBBHomerange <- function(gps, min.fixes, contour.percent=95, output.proj, output.UD=FALSE){
  
  # order the data
  gps <- arrange(gps, AnimalID, acquisitiontime)
  gps$AnimalID <- factor(gps$AnimalID)
  
  # create a table that holds attributes for each animal (to add back to polygons)
  att.table <- gps[match(unique(gps$AnimalID),gps$AnimalID),]
  
  
  # Transform input gps to proj.crs (typically UTM for Oregon)
  gps.sf.utm <- st_transform(gps,crs=output.proj)  # transform to user input or UTM11 WGS84
  
  # convert time to POSIXct
  gps.sf.utm$acquisitiontime <- as.POSIXct(gps.sf.utm$acquisitiontime,tz="UTC")
  
  # convert to spatial points object (sp package) #
  # habitatHR package requires spatial points data frame object
  gps.sp <- as(gps.sf.utm,"Spatial")
  
  print(paste("Kernel function: Brownian bridge"))
  print(paste("Total rows in GPS table for HR calculation: ",nrow(gps.sp)))
  print(paste("Date range for HR calculation: ",(range(gps$acquisitiontime))[1]," to ",(range(gps$acquisitiontime))[2]))
  print(paste("Contour level is set to: ",contour.percent, "%"))
  
  # need to create a ltraj class for the model (note might need to remove duplicates) #
  collar_traj <- as.ltraj(xy = gps.sp@coords, date = gps.sp@data$acquisitiontime, 
                          id = gps.sp@data$AnimalID) 
  
  st <- summary(collar_traj)
  n.collar <- nrow(st)
  
  # # calculate the sig1 smoothing parameter #
  # par(mar=c(1,1,1,1))
  # sig2 <- 30  # sd of relocations as sample of actual position of animal (guess?, collar mfg?)
  # lik <- liker(collar_traj, sig2=sig2, rangesig1=c(0.1,10))
  # # compute average sig1 #
  # tot <- 0
  # for(i in 1:length(lik)){
  #   tot <- tot + lik[[i]][1]$sig1
  # }
  # avg_sig1 <- tot/length(lik)
  
  # based on trials, avg_sig1 is 2.58
  avg_sig1 <- 2.58
  sig2 <- 30  # meters, assumed sd of relocations
  
  
  # BB home range #
  kud.bb <- kernelbb(collar_traj, sig1=avg_sig1 ,sig2=sig2, grid=400, extent=0.5,same4all = TRUE)
  homeranges <- getverticeshr(kud.bb, contour.percent)
  
  # construct output (to keep things simple limiting to either polygon or UD)
  # if (output.UD) output <- kud.bb
  # else output <- homerange_bb
  
  # add some attributes to the polygon data
  newtable <- left_join(homeranges@data,att.table,join_by(id == AnimalID))
  homeranges@data <- newtable[,1:10]
  
  output <- list(homeranges=homeranges,ud=kud.bb)
  
  
  return(output)
}

# this function takes a set of homerange polygons in a projected coord system and returns
# a matrix with the overlap between each pair (assumes they are from the same herd)
calculateHomerangeOverlap <- function(homerange){
  # compute num animals and build matrix
  
    n.animal <- length(homerange$id)
  intersect.mat <- matrix(NA,nrow=n.animal,ncol=n.animal)
  dimnames(intersect.mat) <- list(homerange$id, homerange$id)
  
  homerange <- st_as_sf(homerange)  # convert to sf to use nice intersection function
  
  # loop to compute overlap between each (probably could vectorize this later ?)
  for (i in 1:n.animal){
    for (j in 1:n.animal) {
      
      if (i != j){  # don't overlap with itself!
        ol <- st_intersection(homerange[i,],homerange[j,])
        if (nrow(ol) >0 ) {
          ol$area.overlap_ha <- st_area(ol)/10000
          intersect.mat[i,j] <- ol$area.overlap_ha/ol$area
          intersect.mat[j,i] <- ol$area.overlap_ha/ol$area.1
        }
        else {
          intersect.mat[i,j] <- 0
          intersect.mat[j,i] <- 0
        }
      }
    }
  }
  
  return(intersect.mat)
}

makeGPSMap <- function(data, zcol="AnimalID",colors,alpha=0.8) {
  outmap <- mapview(data, zcol=zcol, legend=TRUE, cex=4,lwd=1, col.regions=colors,alpha=alpha)
  outmap@map
}

# makeHomerangeMap <- function(data, zcol="id",colors,alpha=0.8) {
#   outmap <- mapview(data, zcol=zcol, burst=TRUE,legend=TRUE, cex=4,lwd=1, col.regions=colors,alpha=alpha)
#   outmap@map
# }

makeHomerangeMap <- function(data){
  
  # a large unique color vector for big sets of animalID
    n.animals <- nrow(data@data)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    # handle missing or NA in testing fields so they don't screw up the color mapping
    data@data$CaptureELISAStatus[which(data@data$CaptureELISAStatus==""| is.na(data@data$CaptureELISAStatus),arr.ind=TRUE)] <- "No Record"
    data@data$CapturePCRStatus[which(data@data$CapturePCRStatus==""| is.na(data@data$CapturePCRStatus),arr.ind=TRUE)] <- "No Record"
    
    tcol <- c(brewer.pal(3,'RdYlGn'),"#CCCCCC")
    tcol.shuf <- c(tcol[1:2],tcol[4],tcol[3])
    tcol <- tcol.shuf
    pcr.col.map <- case_when(
      data@data$CapturePCRStatus == "Detected" ~ 1,
      data@data$CapturePCRStatus == "Indeterminate" ~ 2,
      data@data$CapturePCRStatus == "Not detected" ~ 4,
      data@data$CapturePCRStatus == "No Record" ~ 3
    )
    pcr.colors <- tcol[sort(unique(pcr.col.map))]
    #pcr.colors <- tcol[pcr.col.map]
    
    ser.col.map <- case_when(
      data@data$CaptureELISAStatus == "Detected" ~ 1,
      data@data$CaptureELISAStatus == "Indeterminate" ~ 2,
      data@data$CaptureELISAStatus == "Not detected" ~ 4,
      data@data$CaptureELISAStatus == "No Record" ~ 3
    )
    ser.colors <- tcol[sort(unique(ser.col.map))]
    
    n.animals <- length(unique(data@data$id))
    # 
      
      
      if(n.animals > length(col_vector)) id.col <- sample(col_vector, n.animals, replace=TRUE) else id.col <- sample(col_vector, n.animals)
      sex.col <- c("pink","blue")
      celisa.col <- colorRamps::matlab.like(7)
     
      outmap <- mapview(data, zcol="id", burst=FALSE,legend=TRUE, cex=4,lwd=1, col.regions=id.col,alpha=0.8)+
        mapview(data, zcol="Sex", legend=TRUE, cex=4,lwd=1, col.regions=sex.col,alpha=0.8, hide=TRUE)+
        mapview(data, zcol="CapturePCRStatus", legend=TRUE, cex=4,lwd=1, col.regions=pcr.colors,alpha=0.8, hide=TRUE)+
        mapview(data, zcol="CaptureELISAStatus", legend=TRUE, cex=4,lwd=1, col.regions=ser.colors,alpha=0.8, hide=TRUE)+
        mapview(data, zcol="Capture_cELISA", legend=TRUE, cex=4,lwd=1, col.regions=celisa.col,at=c(-15,0,20,40,60,80,100),alpha=0.8, hide=TRUE)
      
    return(outmap@map)
      
}


overlapImagePlot <- function(intersect.mat){
  #   # simple image of the overlap matrix
  # n.col <- ncol(intersect.mat)
  # n.row <- nrow(intersect.mat) 
  # image(1:n.col, 1:n.row, t(intersect.mat), col = brewer.pal(n = 8, name = "YlOrRd"), 
  #       axes = FALSE,xlab='',ylab='', main="Homerange Overlap")
  # par(las=2)
  # axis(1, 1:n.col, colnames(intersect.mat))
  # axis(2, 1:n.row, rownames(intersect.mat))
  # for (x in 1:n.col)
  #   for (y in 1:n.row)
  #     text(x, y, round((intersect.mat)[y,x],2))  
  
  # using the heatmaply package, gives interactive plot and reduced clutter
  heatmaply(signif(intersect.mat,digits=2),Rowv=FALSE,Colv=FALSE,colors= cool_warm,
            column_text_angle=90, label_format_fun=function(...) round(..., digits=2))
}

# function to get nodes and edges from a matrix from DI between animals or overlap of HR between group of animals
overlapClusterDend <- function(data.matrix){
  
  g <- graph_from_adjacency_matrix(data.matrix,mode="directed",weighted=TRUE, diag=FALSE)
  
  #E(g)$width <- E(g)$weight*5
  
  
  # ceb <- cluster_edge_betweenness(g,membership = TRUE)
  # dendPlot(ceb,mode="hclust",main="")
  
  cw <- cluster_walktrap(g,steps=8)
  #cw <- cluster_edge_betweenness(g,directed=TRUE)
  
  V(g)$community <- cw$membership
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  dend <- as.dendrogram(cw)
  nclust <- as.integer(max(cw$membership))
  mycolors1 <- col_vector[1:nclust]
  col_clus <- mycolors1[cw$membership]
  labels_colors(dend) <- col_clus[order.dendrogram(dend)]
  
  plot(dend, main="Cluster Dendrogram of Overlap Network")
  
}

overlapNetworkPlot <- function(data.matrix,gps){
  
  g <- graph_from_adjacency_matrix(data.matrix,mode="directed",weighted=TRUE, diag=FALSE)
  cw <- cluster_walktrap(g,steps=8)
  #cw <- cluster_edge_betweenness(g,directed=TRUE)
  
  V(g)$community <- cw$membership
  
  # qualitative color pallete for groups
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    #mycolors <- sample(col_vector, max(cw$membership))
    mycolors <- col_vector[1:max(cw$membership)]
  
  # first index in gps data frame that matches each animal 
    t.inds <-  match(V(g)$name, gps$AnimalID)
    df <- gps[t.inds,]
    df <- df %>% select(AnimalID, Sex, CapturePCRStatus, CaptureELISAStatus, Capture_cELISA)
    
    V(g)$MoviStatus <- as.factor(df$CaptureELISAStatus)
    V(g)$MoviPCR <- as.factor(df$CapturePCRStatus)
    V(g)$MovicElisa <- df$Capture_cELISA
    
  # colors for test results #
    test.col <- brewer.pal(5,'RdYlGn')
    test.col <- test.col[c(1,3,5)]
    label.col.map <- case_when(
      V(g)$MoviStatus == "Detected" ~ 1,
      V(g)$MoviStatus == "Indeterminate" ~ 2,
      V(g)$MoviStatus == "Not detected" ~ 3
    )
    # set label color to test Elisa results
    V(g)$label.color <- test.col[label.col.map]
  
  # set shape to PCR test 
    shape.map <- case_when(
      V(g)$MoviPCR == "Detected" ~ "star",
      V(g)$MoviPCR == "Indeterminate" ~ "triangle",
      V(g)$MoviPCR == "Not detected" ~ "dot"
    )
    # set label color to test Elisa results
    V(g)$shape <- shape.map
  
  #plot(g, layout=layout, vertex.color=mycolors[V(g)$community],vertex.label=V(g)$name,edge.arrow.size=.2)
 
   E(g)$width <- E(g)$weight*3
  V(g)$label.cex = 0.7
  #V(g)$label.color = "black"
  V(g)$color <- mycolors[V(g)$community]
  V(g)$label=V(g)$name
  visIgraph(g, idToLabel = FALSE) %>% visIgraphLayout(layout = "layout_with_kk") %>%
    visOptions(nodesIdSelection = TRUE, highlightNearest = TRUE)
  
}

getClusterData <- function(data.matrix){
  
  
  # very simple code for this
  g <- graph_from_adjacency_matrix(data.matrix,mode="directed",weighted=TRUE, diag=FALSE)
  
  # we have a variety of cluster options in igraph #
  #ceb <- cluster_edge_betweenness(g,membership = TRUE)
  #co <- cluster_optimal(g)  # not hieracrchical so can't use dendPlot
  cw <- cluster_walktrap(g,steps=8)
  #cw <- cluster_louvain(g)
  #cw <- cluster_edge_betweenness(g,directed=TRUE)
  #dendPlot(ceb,mode="hclust",main="")
  
  #data to return
  return(data.frame(AnimalID=cw$names,Cluster=cw$membership))
  
}

attributeNetworkPlot <- function(data.matrix, display="Elisa", gps){
  
  # create the graph, and cluster
  g <- graph_from_adjacency_matrix(data.matrix,mode="directed",weighted=TRUE, diag=FALSE)
  cw <- cluster_walktrap(g,steps=8)
  #cw <- cluster_edge_betweenness(g,directed=TRUE)
  
  #layout <-layout.fruchterman.reingold(g)
  #layout <- layout_with_kk(g)
  layout <- layout_with_graphopt(g)
  #layout <- layout_with_lgl(g)
  #layout <- layout_nicely(g)
  
  
  V(g)$community <- cw$membership
  
  # join the community to the testing result tables (note this needs to be updated to deal with recaps)
  # df <- data.frame(AnimalID=V(g)$name)
  # p.sub <- pcr %>% filter(AnimalID %in% V(g)$name) %>% select(AnimalID, Sample_ID, SampleDate, MoviPCR)
  # s.sub <- ser %>% filter(AnimalID %in% V(g)$name) %>% select(AnimalID, Sample_ID, SampleDate, Movi_Elisa, Movi_Status)
  # df <- left_join(df,p.sub) %>% left_join(s.sub)
  # 
  
  
  # first index in gps data frame that matches each animal 
  t.inds <-  match(V(g)$name, gps$AnimalID)
  df <- gps[t.inds,]
  df <- df %>% select(AnimalID, Sex, CapturePCRStatus, CaptureELISAStatus, Capture_cELISA)
  
  V(g)$MoviStatus <- as.factor(df$CaptureELISAStatus)
  V(g)$MoviPCR <- as.factor(df$CapturePCRStatus)
  V(g)$MovicElisa <- df$Capture_cELISA
  
  
  # colors for test results #
  tcol <- brewer.pal(3,'RdYlGn')
  
  if (display=="Elisa") {
    
    att.vec <- V(g)$MoviStatus
    legend.title <- "Movi ELISA Status"
    legend <- levels(V(g)$MoviStatus)
    V(g)$group <- V(g)$MoviStatus
    l.title <- "Movi ELISA Status"
  }
  if (display == "PCR") {
    
    att.vec <- V(g)$MoviPCR
    legend.title <- "Movi PCR"
    legend <- levels(V(g)$MoviPCR)
    V(g)$group <- V(g)$MoviPCR
    l.title <- "Movi PCR Status"
  }
  
  # we want same color pattern for test results (having no "Detected" can swap colors between plots)
  col.map <- case_when(
    att.vec == "Detected" ~ 1,
    att.vec == "Indeterminate" ~ 2,
    att.vec == "Not detected" ~ 3
  )
  V(g)$color <- tcol[col.map]
  
  E(g)$width <- E(g)$weight*3
  
  # construct a plot with legend
  #     V(g)$label <- V(g)$name
  # plot(g, layout=layout,edge.arrow.size=.5, 
  #      vertex.label.font=2, vertex.label.color="gray40",
  #      vertex.label.cex=.7, edge.color="gray85")
  # 
  # legend("bottomleft", legend, pch=21,
  #        col="#777777", pt.bg=tcol[sort(unique(col.map))], pt.cex=2, cex=.8, bty="n", ncol=1, title=legend.title)
  # 
  
  visIgraph(g, idToLabel = FALSE) %>% 
    visOptions(nodesIdSelection = TRUE, highlightNearest = TRUE) %>% 
    visIgraphLayout(layout = "layout_nicely") %>%
    visGroups(groupname = "Detected", color = tcol[1]) %>%
    visGroups(groupname = "Indeterminate", color = tcol[2]) %>%
    visGroups(groupname = "Not detected", color = tcol[3]) %>%
    visLegend(enabled=TRUE, useGroups=TRUE, main=paste(l.title), position="left")
  
}

# a static network plot for the download button
overlapNetworkPlotOutput <- function(data.matrix,gps){
  
  g <- graph_from_adjacency_matrix(data.matrix,mode="directed",weighted=TRUE, diag=FALSE)
  cw <- cluster_walktrap(g,steps=8)
  #cw <- cluster_edge_betweenness(g,directed=TRUE)
  
  V(g)$community <- cw$membership
  
  # qualitative color pallete for groups
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #mycolors <- sample(col_vector, max(cw$membership))
  mycolors <- col_vector[1:max(cw$membership)]
  
  # first index in gps data frame that matches each animal 
  t.inds <-  match(V(g)$name, gps$AnimalID)
  df <- gps[t.inds,]
  df <- df %>% select(AnimalID, Sex, CapturePCRStatus, CaptureELISAStatus, Capture_cELISA)
  
  V(g)$MoviStatus <- as.factor(df$CaptureELISAStatus)
  V(g)$MoviPCR <- as.factor(df$CapturePCRStatus)
  V(g)$MovicElisa <- df$Capture_cELISA
  
  # colors for test results #
  test.col <- c("green4","sienna3","red3")
  label.col.map <- case_when(
    V(g)$MoviStatus == "Detected" ~ 1,
    V(g)$MoviStatus == "Indeterminate" ~ 2,
    V(g)$MoviStatus == "Not detected" ~ 3
  )
  # set label color to test Elisa results
  V(g)$label.color <- test.col[label.col.map]
  
  # Define a triangle vertex shape
      mytriangle <- function(coords, v = NULL, params) {
        vertex.color <- params("vertex", "color")
        if (length(vertex.color) != 1 && !is.null(v)) {
          vertex.color <- vertex.color[v]
        }
        vertex.size <- 1 / 200 * params("vertex", "size")
        if (length(vertex.size) != 1 && !is.null(v)) {
          vertex.size <- vertex.size[v]
        }
        
        symbols(
          x = coords[, 1], y = coords[, 2], bg = vertex.color,
          stars = cbind(vertex.size, vertex.size, vertex.size),
          add = TRUE, inches = FALSE
        )
      }
      # clips as a circle
      add_shape("triangle",
                clip = shapes("circle")$clip,
                plot = mytriangle
      )
  
  # set shape to PCR test 
  shape.map <- case_when(
    V(g)$MoviPCR == "Detected" ~ "square",
    V(g)$MoviPCR == "Indeterminate" ~ "triangle",
    V(g)$MoviPCR == "Not detected" ~ "circle"
  )
  # set label color to test Elisa results
  V(g)$shape <- shape.map
  
  
  
  E(g)$width <- E(g)$weight*3
  V(g)$label.cex = 0.7
  #V(g)$label.color = "black"
  V(g)$color <- mycolors[V(g)$community]
  V(g)$label=V(g)$name
  
  
  plot(g, layout=layout_nicely,vertex.labels=V(g)$label)
  legend("bottomleft", legend=c("Detected","Indeterminate","Not detected")  , col = "black" , 
         bty = "n", pch=c(0,2,1) , pt.cex = 1.5, cex = 1.0, text.col="black" , horiz = FALSE,title="PCR Status at capture")
  legend("topleft", legend=c("Detected","Indeterminate","Not detected")  , col = test.col , 
         bty = "n", pch=NA , pt.cex = 1.5, cex = 1.0, text.col=test.col,title.col="black" , horiz = FALSE,title="ELISA Status at capture")
  
}

overlapMatrixPlotDownload <- function(intersect.mat){
    # simple image of the overlap matrix
  n.col <- ncol(intersect.mat)
  n.row <- nrow(intersect.mat)
  
  # a few settings
  if (n.col > 18) ax.cex <- 0.6 else ax.cex <- 1
  
  par(oma=c( 0,0,0,4)) # margin of 4 spaces width at right hand side
  par(mar=c(6,6,3,4)+0.1)
  image(1:n.col, 1:n.row, t(intersect.mat), col = colorspace::diverge_hsv(8),
        axes = FALSE,xlab='',ylab='', main="Homerange Overlap")
  par(las=2)
  axis(1, 1:n.col, colnames(intersect.mat),cex.axis=ax.cex)
  axis(2, 1:n.row, rownames(intersect.mat),cex.axis=ax.cex)
  if (n.col < 21) { # only draw text labels on small graphs
    for (x in 1:n.col)
    for (y in 1:n.row)
      text(x, y, round((intersect.mat)[y,x],1),col='grey50')
  }

  par(oma=c( 0,0,0,1))# reset margin to be much smaller.
  image.plot( legend.only=TRUE,zlim=c(0,1), col=colorspace::diverge_hsv(8)) 
}


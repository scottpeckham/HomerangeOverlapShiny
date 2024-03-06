# Helper file for Homerange overlap v2 Shiny app
# this script houses some functions to compute homeranges from GPS tables from our DB#
#   assumes lat/lon are in WGS84 (i.e. vectronics data format)

# required packages #
library(sf)
library(RColorBrewer)
library(mapview)
library(dplyr)
library(igraph)
library(visNetwork)
library(colorRamps)
library(fields)
library(colorspace)
# required packages #



# a static network plot for the download button
overlapNetworkPlotOutput <- function(community, display="Both", gps){
  
  # create the graph, and cluster
  g <- community$graph
  cw <- community$cluster
  
  V(g)$community <- cw$membership
  
  # qualitative color pallete for groups
    # lots of colors for IDs (if needed)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    # better paired pallette for clusters
    clus.colors <- brewer.pal(10,'Paired') 
  
    #mycolors <- sample(col_vector, max(cw$membership))
    mycolors <- clus.colors[1:max(cw$membership)]
  
  # first index in gps data frame that matches each animal 
  t.inds <-  match(V(g)$name, gps$AnimalID)
  df <- gps[t.inds,]
  df <- df %>% select(AnimalID, Sex, CapturePCRStatus, CaptureELISAStatus, Capture_cELISA)
  
  V(g)$MoviStatus <- as.factor(df$CaptureELISAStatus)
  V(g)$MoviPCR <- as.factor(df$CapturePCRStatus)
  V(g)$MovicElisa <- df$Capture_cELISA
  
  # colors for test results #
  test.col <- c("green4","sienna3","red3","grey70")
  label.col.map <- case_when(
    V(g)$MoviStatus == "Detected" ~ 1,
    V(g)$MoviStatus == "Indeterminate" ~ 2,
    V(g)$MoviStatus == "Not detected" ~ 3,
    .default = 4
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
    V(g)$MoviPCR == "Not detected" ~ "circle",
    .default = "circle"
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

getClusterData <- function(community){
  
  
  # very simple code for this, just strip from community list object and make a df
  cw <- community$cluster
  
  #data to return
  return(data.frame(AnimalID=cw$names,Cluster=cw$membership))
  
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

ClusterMapview <- function(hr,idcolors,map.clus.colors){
  m <- mapview(hr,zcol="AnimalID",col.regions=idcolors,alpha.regions=0.75)+ 
    mapview(hr, zcol="Cluster", legend=TRUE, cex=2,lwd=1, col.regions=map.clus.colors,alpha.regions=0.8)+ 
    mapview(hr,zcol="cluster.mean.cELISA")+ 
    mapview(hr,zcol="cluster.ELISA.prev")+
    mapview(hr,zcol="cluster.PCR.prev")
}


# this script houses the workhorse functions commonly used to compute homeranges from GPS tables from our DB#
#   and in many other codes that deal with home ranges, clustering, and overlap
#   assumes lat/lon are in WGS84  EPSG 4326

# required packages #
library(sf)
library(adehabitatHR)
library(adehabitatLT)
library(RColorBrewer)
library(mapview)
library(odbc)
library(RSQLite)
library(sp)
library(dplyr)
library(igraph)
library(heatmaply)
library(dendextend)
library(visNetwork)
library(sfdep)

# required packages #

# Calculate a biological year vector from a date vector
  getBioYear <- function(indate){
    year <- as.numeric(strftime(indate, format="%Y",tz="UTC"))
    month <- as.numeric(strftime(indate, format = "%m",tz="UTC"))
    if(month > 4) bioyear <- year else bioyear <- year-1
    return(bioyear)
  }

# Add biological year to a data frame (calls getBioYear above)
  addBioYear <- function(gps.data){
    # create a capture year column #
    gps.data$BioYear <-unlist(lapply(gps.data$acquisitiontime, getBioYear))
    return(gps.data)
  }

# For inside the BB home range function
  foo <- function(dt){
    return(dt> (7*3600*24))
    }
# Calculate Animal home range using the standard href kernel 
#
# Inputs:
#       gps - sf object containing points to compute over, with AnimalID the first column, spatial ref EPSG 4326
#       min.fixes - minimum number of locations an animal can have, skip otherwise
#       contour.percent - the percentage used in adehabitat 'getverticesHR' function on the UD, default to 95
#       output.proj - integer denoting appropriate EPSG CRS for the data, typically a UTM, State Plane, or other projected CRS
#       output.UD - logical, should the UD be output instead of home range polygons
#       grid - integer, size of grid used (see adehabitatHR) default to 800
#       extent - numeric, extent of region used in 'kernelUD' (see adehabitatHR), default to 0.5 but this may need to be 
#                 ~1 for 95% home ranges with large dispersion, most common cause for 'getverticesHR' function failing
# Output: 
#       Either a sf object of home ranges, or a named list with the homeranges and estUD object depending on choice above
  
calculateHomerange <- function(gps, min.fixes, output.proj, contour.percent=95, grid=800, extent=1.0, output.UD=FALSE){
  
  # check data type
  if(!is(gps, "sf"))
    stop("first argument has to be a spatial object ")
  
  # check number of locations per animal, adehabitat requires >5/per
  loc.per.animal <- gps %>% st_drop_geometry() %>% group_by(AnimalID) %>%
    summarize(Fix_count =  n())
  drops <- which(loc.per.animal$Fix_count < min.fixes,arr.ind=TRUE)
  excludes <- loc.per.animal$AnimalID[drops]
  if (length(drops)>0) {
    gps <- gps[!gps$AnimalID %in% excludes,]
    print(paste("Excluded",length(drops),"animals with fix count <",min.fixes))
  }
  
  # create a table that holds attributes for each animal (to add back to polygons)
  att.table <- gps[match(unique(gps$AnimalID),gps$AnimalID),] %>% as.data.frame()
  sz <- gps %>% as.data.frame() %>% group_by(AnimalID) %>% summarize(nLocations=n(),startDate=min(acquisitiontime),endDate=max(acquisitiontime))
  att.table <- left_join(att.table,sz,by="AnimalID")
  #att.table <- att.table[,c(1:8,10,12:14)] #strip off some useless columns
  
  # remove any possible duplicates
  gps <- gps %>% distinct()
  
  # ensure AnimalID is a factor
  gps$AnimalID <- factor(gps$AnimalID)
  gps.sf <- subset(gps,select= AnimalID) # Keep just the "AnimalID" column of spatial object (again a restraint of the package) 
  
  # Transform input gps to proj.crs (typically UTM for Oregon)
  gps.sf <- st_transform(gps.sf,crs=output.proj)  # transform to user input
  
  # convert from sf to spatial points data frame (reqd by adehabitatHR)
  gps.collarIDT <- as(gps.sf, "Spatial")
  
  print(paste("Kernel function: Bivariate normal"))
  print(paste("Total number of rows in GPS table for HR calculation: ",nrow(gps)))
  print(paste("Date range for HR calculation: ",(range(gps$acquisitiontime))[1]," to ",(range(gps$acquisitiontime))[2]))
  print(paste("Contour level is set to: ",contour.percent, "%"))
  
  
  # compute UD and polygon HR
  kud <- kernelUD(gps.collarIDT, grid=grid,same4all = FALSE, h="href",extent=extent)
  homeranges <- getverticeshr(kud, percent=contour.percent)
  
  # add some attributes to the polygon data
  newtable <- left_join(homeranges@data,att.table,join_by(id == AnimalID))
  homeranges@data <- newtable
  
  homeranges.out <- st_as_sf(homeranges)   # convert to sf, for some reason the BB homeranges CRS are screwed up by adehabitat, 
  st_crs(homeranges.out) <- output.proj
  
  # Restore column name
  names(homeranges.out)[1] <- "AnimalID"
  
  nanimals <- nrow(homeranges.out)
  print(paste("Number of animals: ",nanimals))
  
  # construct output (to keep things simple limiting to either polygon or UD)
  if (output.UD) output <- list(homeranges=homeranges.out,ud=kud)
  else output <- homeranges.out
  
  
  return(output)
  
}

# Calculate Animal home range using the Brownian Bridge kernel 
#
# Inputs:
#       gps - sf object containing points to compute over, with AnimalID the first column, spatial ref EPSG 4326
#       min.fixes - minimum number of locations an animal can have, skip otherwise
#       contour.percent - the percentage used in adehabitat 'getverticesHR' function on the UD, default to 95
#       output.proj - integer denoting appropriate EPSG CRS for the data, typically a UTM, State Plane, or other projected CRS
#       output.UD - logical, should the UD be output instead of home range polygons
#       grid - integer, size of grid used (see adehabitatHR) default to 800
#       extent - numeric, extent of region used in 'kernelUD' (see adehabitatHR), default to 0.5 but this may need to be 
#                 ~1 for 95% home ranges with large dispersion, most common cause for 'getverticesHR' function failing
# Output: 
#       Either a sf object of home ranges, or estUD object depening on choice above
calculateBBHomerange <- function(gps, min.fixes, contour.percent=95, output.proj, output.UD=FALSE, grid=800, extent=1.0){
  
  # check data type
  if(!is(gps, "sf"))
    stop("first argument has to be a spatial object ")
  
  # check number of locations per animal, adehabitat requires >5/per
  #   and exclude any with < min.fixes
  loc.per.animal <- gps %>% st_drop_geometry() %>% group_by(AnimalID) %>%
    summarize(Fix_count =  n())
  drops <- which(loc.per.animal$Fix_count < min.fixes,arr.ind=TRUE)
  excludes <- loc.per.animal$AnimalID[drops]
  if (length(drops)>0) {
    gps <- gps[!gps$AnimalID %in% excludes,]
    print(paste("Excluded",length(drops),"animals with fix count <",min.fixes))
  }
  
  # order the data
  gps <- arrange(gps, AnimalID, acquisitiontime)
  gps$AnimalID <- factor(gps$AnimalID)
  
  # remove any possible duplicates
  gps <- gps %>% distinct()
  
  # create a table that holds attributes for each animal (to add back to polygons)
  att.table <- gps[match(unique(gps$AnimalID),gps$AnimalID),] 
  sz <- gps %>% as.data.frame() %>% group_by(AnimalID) %>% summarize(nLocations=n(),startDate=min(acquisitiontime),endDate=max(acquisitiontime))
  att.table <- left_join(att.table,sz,by="AnimalID")
  #att.table <- att.table[,c(1:8,10,12:14)] #strip off some useless columns
  
  # Transform input gps to proj.crs (typically UTM for Oregon, but must be projected coord system)
  gps.sf.utm <- st_transform(gps,crs=output.proj)  # transform to user input 
  
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
  
  # Cut the trajectory to deal with large "gaps" in the data caused by ? or when an animal switches collars 
    # function to return dt > 30 days
    
    collar_traj <- cutltraj(collar_traj, "foo(dt)", nextr = TRUE)
    
  st <- summary(collar_traj)
  n.collar <- length(unique(st$id))
  print(paste("Number of individuals in the data set: ",n.collar))
  
  # # calculate the sig1 smoothing parameter, brownian motion variance parameter is s1^2 #
    # par(mar=c(1,1,1,1))
    # sig2 <- 50  # sd of relocations as sample of actual position of animal (guess?, collar mfg?)
    # lik <- liker(collar_traj, sig2=sig2, rangesig1=c(0.1,20))
    # # compute average sig1 #
    # tot <- 0
    # for(i in 1:length(lik)){
    #   tot <- tot + lik[[i]][1]$sig1
    # }
    # avg_sig1 <- tot/length(lik)
  
  # based on trials, avg_sig1 is 2.58 for LM, 2.2 for LS
  avg_sig1 <- 2.58
  sig2 <- 30  # meters, assumed sd of relocations
  
  
  # BB home range #
  kud.bb <- kernelbb(collar_traj, sig1=avg_sig1 ,sig2=sig2, grid=grid, extent=extent,same4all = FALSE)
  homerange_bb <- getverticeshr(kud.bb, contour.percent)
  
  # add some attributes to the polygon data
  newtable <- left_join(homerange_bb@data,att.table,join_by(id == AnimalID))
  homerange_bb@data <- newtable
  
  homerange_bb.out <- st_as_sf(homerange_bb)   # for some reason the BB homeranges CRS are screwed up by adehabitat
  st_crs(homerange_bb.out) <- output.proj
  
  # Restore column name
  names(homerange_bb.out)[1] <- "AnimalID"
  
  # construct output (if UD=TRUE we'll output a list ojecte, if FALSE (default) output is sf polygons)
  if (output.UD) output <- list(homeranges=homerange_bb.out, ud=kud.bb)
  else output <- homerange_bb.out
  
  
  return(output)
}


# this function takes a set of homerange polygons in a projected coord system and returns
# a matrix with the overlap between each pair (assumes they are from the same herd)
# **** assumes animal ID is column 1 and named 'AnimalID'
calculateHomerangeOverlap <- function(homerange){
  # compute num animals and build matrix
  
  # check data type
  if(!is(homerange, "sf"))
    stop("first argument has to be a spatial object of type 'sf' ")
  
  n.animal <- length(homerange$AnimalID)
  intersect.mat <- matrix(NA,nrow=n.animal,ncol=n.animal)
  dimnames(intersect.mat) <- list(homerange$AnimalID, homerange$AnimalID)
  
  
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

# Same function as above but does overlap of Clusters (did this to be quick and easy, should have mod'd the above but this is descriptive)
calculateClusterOverlap <- function(cluster.sf){
  # compute num animals and build matrix
  
  # check data type
  if(!is(cluster.sf, "sf"))
    stop("first argument has to be a spatial object of type 'sf' ")
  
  n.animal <- length(cluster.sf$Cluster)
  intersect.mat <- matrix(NA,nrow=n.animal,ncol=n.animal)
  dimnames(intersect.mat) <- list(cluster.sf$Cluster, cluster.sf$Cluster)
  
  
  # loop to compute overlap between each (probably could vectorize this later ?)
  for (i in 1:n.animal){
    for (j in 1:n.animal) {
      
      if (i != j){  # don't overlap with itself!
        ol <- st_intersection(cluster.sf[i,],cluster.sf[j,])
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

# retained this function for now, but typically used in Shiny app programming
makeGPSMap <- function(data, zcol="AnimalID",colors,alpha=0.8) {
  outmap <- mapview(data, zcol=zcol, legend=TRUE, cex=4,lwd=1, col.regions=colors,alpha=alpha)
  outmap@map
}

# makeHomerangeMap <- function(data, zcol="id",colors,alpha=0.8) {
#   outmap <- mapview(data, zcol=zcol, burst=TRUE,legend=TRUE, cex=4,lwd=1, col.regions=colors,alpha=alpha)
#   outmap@map
# }
makeHomerangeMap <- function(data, colors){
  
  # handle missing or NA in testing fields so they don't screw up the color mapping
  data$CaptureELISAStatus[which(data$CaptureELISAStatus==""| is.na(data$CaptureELISAStatus),arr.ind=TRUE)] <- "No Record"
  data$CapturePCRStatus[which(data$CapturePCRStatus==""| is.na(data$CapturePCRStatus),arr.ind=TRUE)] <- "No Record"
  
  # we want same color pattern for test results (having no "Detected" can swap colors between plots)
  # colors for test results #
  tcol <- c(brewer.pal(3,'RdYlGn'),"#CCCCCC")
  tcol.shuf <- c(tcol[1:2],tcol[4],tcol[3])
  tcol <- tcol.shuf
  pcr.col.map <- case_when(
    data$CapturePCRStatus == "Detected" ~ 1,
    data$CapturePCRStatus == "Indeterminate" ~ 2,
    data$CapturePCRStatus == "Not detected" ~ 3,
    data$CapturePCRStatus == "No Record" ~ 4
  )
  pcr.colors <- tcol[sort(unique(pcr.col.map))]
  #pcr.colors <- tcol[pcr.col.map]
  
  ser.col.map <- case_when(
    data$CaptureELISAStatus == "Detected" ~ 1,
    data$CaptureELISAStatus == "Indeterminate" ~ 2,
    data$CaptureELISAStatus == "Not detected" ~ 4,
    data$CaptureELISAStatus == "No Record" ~ 3
  )
  ser.colors <- tcol[sort(unique(ser.col.map))]
  
  
  id.col <- colors
  sex.col <- c("pink","blue")
  celisa.col <- colorRamps::matlab.like(7)
  
  outmap <- mapview(data, zcol="AnimalID", burst=FALSE,legend=TRUE, cex=4,lwd=1, col.regions=id.col,alpha=0.8)+
    mapview(data, zcol="Sex", legend=TRUE, cex=4,lwd=1, col.regions=sex.col,alpha=0.8)+
    mapview(data, zcol="CapturePCRStatus", legend=TRUE, cex=4,lwd=1, col.regions=pcr.colors,alpha=0.8)+
    mapview(data, zcol="CaptureELISAStatus", legend=TRUE, cex=4,lwd=1, col.regions=ser.colors,alpha=0.8)+
    mapview(data, zcol="Capture_cELISA", legend=TRUE, cex=4,lwd=1, col.regions=celisa.col,at=c(-15,0,20,40,60,80,100),alpha=0.8)
  
  return(outmap@map)
  
}

makeOverviewMap <- function(data,pts,colors,alpha=0.8){
  # handle missing or NA in testing fields so they don't screw up the color mapping
  data$CaptureELISAStatus[which(data$CaptureELISAStatus==""| is.na(data$CaptureELISAStatus),arr.ind=TRUE)] <- "No Record"
  data$CapturePCRStatus[which(data$CapturePCRStatus==""| is.na(data$CapturePCRStatus),arr.ind=TRUE)] <- "No Record"
  
  # we want same color pattern for test results (having no "Detected" can swap colors between plots)
  # colors for test results #
  tcol <- c(brewer.pal(3,'RdYlGn'),"#CCCCCC")
  tcol.shuf <- c(tcol[1:2],tcol[4],tcol[3])
  tcol <- tcol.shuf
  pcr.col.map <- case_when(
    data$CapturePCRStatus == "Detected" ~ 1,
    data$CapturePCRStatus == "Indeterminate" ~ 2,
    data$CapturePCRStatus == "Not detected" ~ 3,
    data$CapturePCRStatus == "No Record" ~ 4
  )
  pcr.colors <- tcol[sort(unique(pcr.col.map))]
  #pcr.colors <- tcol[pcr.col.map]
  
  ser.col.map <- case_when(
    data$CaptureELISAStatus == "Detected" ~ 1,
    data$CaptureELISAStatus == "Indeterminate" ~ 2,
    data$CaptureELISAStatus == "Not detected" ~ 4,
    data$CaptureELISAStatus == "No Record" ~ 3
  )
  ser.colors <- tcol[sort(unique(ser.col.map))]
  
  
  id.col <- colors
  sex.col <- c("pink","blue")
  celisa.col <- colorRamps::matlab.like(7)
  
  outmap <- mapview(data, zcol="AnimalID", burst=FALSE,legend=TRUE, cex=4,lwd=1, col.regions=id.col,alpha=0.8)+
    mapview(data, zcol="Sex", legend=TRUE, cex=4,lwd=1, col.regions=sex.col,alpha=0.8)+
    mapview(data, zcol="CapturePCRStatus", legend=TRUE, cex=4,lwd=1, col.regions=pcr.colors,alpha=0.8)+
    mapview(data, zcol="CaptureELISAStatus", legend=TRUE, cex=4,lwd=1, col.regions=ser.colors,alpha=0.8)+
    mapview(data, zcol="Capture_cELISA", legend=TRUE, cex=4,lwd=1, col.regions=celisa.col,at=c(-15,0,20,40,60,80,100),alpha=0.8)+
    mapview(pts,zcol="AnimalID",col.regions=id.col, cex=3)
  
  return(outmap@map)
  
}
# Use heatmaply package to make quick matrix plot, or uncomment lower code for using 'image'
overlapImagePlot <- function(intersect.mat){
  
  
  # using the heatmaply package, gives interactive plot and reduced clutter
  heatmaply(signif(intersect.mat,digits=2),Rowv=FALSE,Colv=FALSE,colors= cool_warm,
            column_text_angle=90, label_format_fun=function(...) round(..., digits=2))
  
  # # simple image of the overlap matrix (replaced with heatmaply package above)
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
}

# Plot and return a dendrogram from a community list object created by getClusterCommunity
overlapClusterDend <- function(c.list){
  
  # create the graph, and cluster
  g <- c.list$graph
  cw <- c.list$cluster
  
  # assign cluster membership to the community structure
  V(g)$community <- cw$membership
  
  # qual pallette for large numbers of animals (if needed)
    #qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    #col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # better paired pallette for clusters
    clus.colors <- brewer.pal(12,'Paired') 
  
  #mycolors <- sample(col_vector, max(cw$membership))
  #mycolors <- clus.colors[1:max(cw$membership)]
  
  dend <- as.dendrogram(cw)
  nclust <- as.integer(max(cw$membership))
  #mycolors1 <- col_vector[1:nclust]  # using the qual pallette always gave terribly light yellow and some others
  mycolors1 <- clus.colors[1:nclust]
  col_clus <- mycolors1[cw$membership]
  labels_colors(dend) <- col_clus[order.dendrogram(dend)]
  
  if ( length(cw$membership)> 50) dend <- set(dend, "labels_cex", 0.7) # let's decrease the size of the labels when N is large
  
  plot(dend, main="Cluster Dendrogram of Overlap Network")
  
  return(dend) # for use elsewhere
  
}

# Plot community network using igraph built in plotting routine
overlapNetworkPlot <- function(c.list){
  
  # create the graph, and cluster
  g <- c.list$graph
  cw <- c.list$cluster
  
  layout <- layout_nicely(g)
  
  
  V(g)$community <- cw$membership
  
  # qual pallette for large N
    # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    # mycolors <- col_vector[1:max(cw$membership)]
  
  # better paired pallette for clusters
    clus.colors <- brewer.pal(12,'Paired') 
    mycolors <- clus.colors[1:max(cw$membership)]
  
  #plot(g, layout=layout, vertex.color=mycolors[V(g)$community],vertex.label=V(g)$name,edge.arrow.size=.2)
  
  V(g)$label.cex = 0.7
  V(g)$label.color = "black"
  V(g)$color <- mycolors[V(g)$community]
  V(g)$label=V(g)$name
  visIgraph(g, idToLabel = FALSE) %>% visIgraphLayout(layout = "layout_nicely") %>%
    visOptions(nodesIdSelection = TRUE, highlightNearest = TRUE)
  
  #plot(cw, g, layout=layout, vertex.size=5,  vertex.label=V(g)$name, edge.arrow.size=.2)
  #plot(cw,g)
  
}

# Compute our igraph object and detect clusters
# function to get nodes and edges from a matrix from DI between animals or overlap of HR between group of animals
#   current available methods for clustering community are "walktrap", "louvain", "edge_betweenness" see igprah page for detail
#   only walktrap and edge_b are avaialbe for directed adjacency matrices, louvain is not
getClusterCommunity <- function(data.matrix,method="walktrap", mode="directed", steps=4){
  
  
  g <- graph_from_adjacency_matrix(data.matrix,mode=mode,weighted=TRUE, diag=FALSE)
  
  if (method=="walktrap"){
    cw <- cluster_walktrap(g,steps=steps)
  }
  if (method=="louvain"){
    cw <- cluster_louvain(g)
  }
  if (method=='edge_betweenness'){
    cw <- cluster_edge_betweenness(g,membership=TRUE)
  }
  
  # assign cluster membership to the community structure
  V(g)$community <- cw$membership
  
  output <- list(graph=g, cluster=cw)
  #data to return
  return(output)
  
}

getClusterDataFrame <- function(cluster){
  
  
  #data to return
  return(data.frame(AnimalID=cluster$names,Cluster=cluster$membership))
  
}

# Make a attributed network plot #
# input is a list object from getClusterCommunity with names 'graph' and 'cluster'
attributeNetworkPlot <- function(c.list, display="Elisa", gps){
  
  # create the graph, and cluster
  g <- c.list$graph
  cw <- c.list$cluster
  
  V(g)$community <- cw$membership
  
  # join the community to the testing result tables (note this needs to be updated to deal with recaps)
  # df <- data.frame(AnimalID=V(g)$name)
  # p.sub <- pcr %>% filter(AnimalID %in% V(g)$name) %>% select(AnimalID, Sample_ID, SampleDate, MoviPCR)
  # s.sub <- ser %>% filter(AnimalID %in% V(g)$name) %>% select(AnimalID, Sample_ID, SampleDate, Movi_Elisa, Movi_Status)
  # df <- left_join(df,p.sub) %>% left_join(s.sub)
  # 
  
  # qualitative color pallete for groups or large N
    # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    # #mycolors <- sample(col_vector, max(cw$membership))
    # mycolors <- col_vector[1:max(cw$membership)]
  
  # better paired pallette for clusters
    clus.colors <- brewer.pal(12,'Paired') 
    mycolors <- clus.colors[1:max(cw$membership)]
  
  # first index in gps data frame that matches each animal (note this can also be sf homerange poly with attributes)
  t.inds <-  match(V(g)$name, gps$AnimalID)
  df <- gps[t.inds,]
  df <- df %>% select(AnimalID, Sex, CapturePCRStatus, CaptureELISAStatus, Capture_cELISA)
  
  V(g)$MoviStatus <- as.factor(df$CaptureELISAStatus)
  V(g)$MoviPCR <- as.factor(df$CapturePCRStatus)
  V(g)$MovicElisa <- df$Capture_cELISA
  
  
  # colors for test results #
  tcol <- brewer.pal(5,'RdYlGn')
  tcol <- brewer.pal(5,'RdBu')
  tcol <- tcol[c(1,3,5)]
  
  if (display=="Elisa") {
    
    att.vec <- V(g)$MoviStatus
    legend.title <- "Movi ELISA Status"
    legend <- levels(V(g)$MoviStatus)
    V(g)$group <- V(g)$MoviStatus
    l.title <- "Movi ELISA Status"
    
    # we want same color pattern for test results (having no "Detected" can swap colors between plots)
    col.map <- case_when(
      att.vec == "Detected" ~ 1,
      att.vec == "Indeterminate" ~ 2,
      att.vec == "Not detected" ~ 3,
      .default = NA
    )
    V(g)$color <- tcol[col.map]
  }
  if (display == "PCR") {
    
    att.vec <- V(g)$MoviPCR
    legend.title <- "Movi PCR"
    legend <- levels(V(g)$MoviPCR)
    V(g)$group <- V(g)$MoviPCR
    l.title <- "Movi PCR Status"
    
    # we want same color pattern for test results (having no "Detected" can swap colors between plots)
    col.map <- case_when(
      att.vec == "Detected" ~ 1,
      att.vec == "Indeterminate" ~ 2,
      att.vec == "Not detected" ~ 3,
      .default = NA
    )
    V(g)$color <- tcol[col.map]
  }
  
  if (display == "Both"){
    # colors for test results #
    test.col <- tcol
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
    # set shape to PCR results
    V(g)$shape <- shape.map
    V(g)$group <- V(g)$community
    V(g)$color <- mycolors[V(g)$community]
    
    l.title <- "Community Graph with Movi PCR and ELISA Status"
  } 
  
  
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
  
  if (display == "Both"){
    visIgraph(g, idToLabel = FALSE) %>% 
      visOptions(nodesIdSelection = TRUE, highlightNearest = TRUE) %>% 
      visIgraphLayout(layout = "layout_nicely") #%>%
    #visLegend(enabled=TRUE, useGroups=TRUE, main=paste(l.title), position="left")
    # visGroups(groupname = "Detected", color = tcol[1]) %>%
    # visGroups(groupname = "Indeterminate", color = tcol[2]) %>%
    # visGroups(groupname = "Not detected", color = tcol[3]) %>%
    # visLegend(enabled=TRUE, useGroups=TRUE, main=paste(l.title), position="left")
    
  } else {
    visIgraph(g, idToLabel = FALSE) %>% 
      visOptions(nodesIdSelection = TRUE, highlightNearest = TRUE) %>% 
      visIgraphLayout(layout = "layout_nicely") %>%
      visGroups(groupname = "Detected", color = tcol[1]) %>%
      visGroups(groupname = "Indeterminate", color = tcol[2]) %>%
      visGroups(groupname = "Not detected", color = tcol[3]) %>%
      visLegend(enabled=TRUE, useGroups=TRUE, main=paste(l.title), position="left")
  }
}

#  #
removeMissingGPS <- function(gps){
  # check for missing lat/lon and drop #
  missing.loc <- which(is.na(gps$Latitude),arr.ind=TRUE)
  if (length(missing.loc) > 0) gps <- gps[-(missing.loc),]
  return(gps)
}

# add cluster stats to data frame in sf object
addClusterStats <- function(h.poly) {
  
  # store status as Factor so we can get 0 counts
  ef <- factor(h.poly$CaptureELISAStatus,levels=c("Detected","Not detected", "Indeterminate"))
  pf <- factor(h.poly$CapturePCRStatus,levels=c("Detected","Not detected", "Indeterminate"))
  
  df <- data.frame(Cluster=h.poly$Cluster,ef=ef,pf=pf)
  counts <- df %>% group_by(Cluster) %>% summarise(n=n())
  
  cElisa.count <- df %>% group_by(Cluster) %>% count(ef,.drop=FALSE)
  cPCR.count <- df %>% group_by(Cluster) %>% count(pf,.drop=FALSE)
  
  clusters <- sort(unique(h.poly$Cluster))
  nclusters <- length(clusters)
  cluster.PCR.prev <- c(rep(0,nclusters))
  cluster.ELISA.prev <- c(rep(0,nclusters))
  cluster.nMembers <- c(rep(0,nclusters))
  for (i in 1:nclusters){
    nde <- cElisa.count %>% filter(Cluster==clusters[i] & ef=="Detected") %>% select(n)
    ndp <- cPCR.count %>% filter(Cluster==clusters[i] & pf=="Detected") %>% select(n)
    tot <- counts %>% filter(Cluster==clusters[i]) %>% select(n)
    cluster.PCR.prev[clusters[i]] <- ndp$n/tot$n
    cluster.ELISA.prev[clusters[i]] <- nde$n/tot$n
    cluster.nMembers[clusters[i]] <- tot$n
  }
  df2 <- data.frame(cluster=clusters, cluster.PCR.prev,cluster.ELISA.prev, nMembers=cluster.nMembers)
  
  h.poly$cluster.PCR.prev <- c(rep(0,nrow(h.poly)))
  h.poly$cluster.ELISA.prev <- c(rep(0,nrow(h.poly)))
  h.poly$cluster.mean.cELISA <- c(rep(0,nrow(h.poly)))
  h.poly$cluster.nMembers <- c(rep(0,nrow(h.poly)))
  
  cmean <- h.poly %>% st_drop_geometry() %>% group_by(Cluster) %>% summarise(clus_mean_cElisa=mean(Capture_cELISA,na.rm=TRUE))
  
  for (i in 1:nclusters){
    clus <- clusters[i]
    ind1 <- which(h.poly$Cluster %in% clus, arr.ind=TRUE)
    ind2 <- which(df2$cluster %in% clus, arr.ind=TRUE)
    ind3 <- which(cmean$Cluster %in% clus, arr.ind=TRUE)
    pcr.value <- df2$cluster.PCR.prev[ind2]
    el.value <- df2$cluster.ELISA.prev[ind2]
    m.value <- cmean$clus_mean_cElisa[ind3]
    n.value <- df2$nMembers[ind2]
    h.poly$cluster.PCR.prev[ind1] <- c(rep(pcr.value,length(ind1)))
    h.poly$cluster.ELISA.prev[ind1] <- c(rep(el.value,length(ind1)))
    h.poly$cluster.mean.cELISA[ind1] <- c(rep(m.value,length(ind1)))
    h.poly$cluster.nMembers[ind1] <- c(rep(n.value,length(ind1)))
    
  }
  
  return(h.poly)
}

# add cluster stats to data frame in sf object
addClusterStats2 <- function(h.poly, entryYear.filter=NULL) {
  
  dat <- h.poly %>% st_drop_geometry() # get data table
  
   # if we only want to use a subste based on EntryBioYear to compute stats #
    if (length(entryYear.filter)>0){
      dat <- dat %>% filter(EntryBioYear==entryYear.filter)
    }
  # store status as Factor so we can get 0 counts
  ef <- factor(dat$CaptureELISAStatus,levels=c("Detected","Not detected", "Indeterminate"))
  pf <- factor(dat$CapturePCRStatus,levels=c("Detected","Not detected", "Indeterminate"))
  
  df <- data.frame(Cluster=dat$Cluster,ef=ef,pf=pf)
  counts <- df %>% group_by(Cluster) %>% summarise(n=n())
  
  cElisa.count <- df %>% group_by(Cluster) %>% count(ef,.drop=FALSE)
  cPCR.count <- df %>% group_by(Cluster) %>% count(pf,.drop=FALSE)
  
  clusters <- sort(unique(dat$Cluster))
  nclusters <- length(clusters)
  cluster.PCR.prev <- c(rep(0,nclusters))
  cluster.ELISA.prev <- c(rep(0,nclusters))
  cluster.nMembers <- c(rep(0,nclusters))
  for (i in 1:nclusters){
    nde <- cElisa.count %>% filter(Cluster==clusters[i] & ef=="Detected") %>% select(n)
    ndp <- cPCR.count %>% filter(Cluster==clusters[i] & pf=="Detected") %>% select(n)
    tot <- counts %>% filter(Cluster==clusters[i]) %>% select(n)
    cluster.PCR.prev[clusters[i]] <- ndp$n/tot$n
    cluster.ELISA.prev[clusters[i]] <- nde$n/tot$n
    cluster.nMembers[clusters[i]] <- tot$n
  }
  df2 <- data.frame(cluster=clusters, cluster.PCR.prev,cluster.ELISA.prev, nMembers=cluster.nMembers)
  
  h.poly$cluster.PCR.prev <- c(rep(0,nrow(h.poly)))
  h.poly$cluster.ELISA.prev <- c(rep(0,nrow(h.poly)))
  h.poly$cluster.mean.cELISA <- c(rep(0,nrow(h.poly)))
  h.poly$cluster.nMembers <- c(rep(0,nrow(h.poly)))
  
  cmean <- h.poly %>% st_drop_geometry() %>% group_by(Cluster) %>% summarise(clus_mean_cElisa=mean(Capture_cELISA,na.rm=TRUE))
  
  for (i in 1:nclusters){
    clus <- clusters[i]
    ind1 <- which(h.poly$Cluster %in% clus, arr.ind=TRUE)
    ind2 <- which(df2$cluster %in% clus, arr.ind=TRUE)
    ind3 <- which(cmean$Cluster %in% clus, arr.ind=TRUE)
    pcr.value <- df2$cluster.PCR.prev[ind2]
    el.value <- df2$cluster.ELISA.prev[ind2]
    m.value <- cmean$clus_mean_cElisa[ind3]
    n.value <- df2$nMembers[ind2]
    h.poly$cluster.PCR.prev[ind1] <- c(rep(pcr.value,length(ind1)))
    h.poly$cluster.ELISA.prev[ind1] <- c(rep(el.value,length(ind1)))
    h.poly$cluster.mean.cELISA[ind1] <- c(rep(m.value,length(ind1)))
    h.poly$cluster.nMembers[ind1] <- c(rep(n.value,length(ind1)))
    
  }
  
  return(h.poly)
}

# Draw trajectory and points, note that I haven't solved how to get the colors to match up correctly w/out summarizing point data and 
#   losing the attributes because I needed to create a linestring using "group_by" - and had to create a multipoint to match using "group_by" so point attributes are lost
makeLinePointMap <- function(sf.dat){
  
  trajectory <- sf.dat %>%
   group_by(AnimalID) %>%
   dplyr::summarize() %>%  
    st_cast("LINESTRING") %>% arrange(AnimalID)
 #points <- sf.dat %>% arrange(AnimalID)

    #trajectory <- trajectory[1:15,]
    
 points <- sf.dat %>%
   group_by(AnimalID) %>%
   dplyr::summarize() %>% arrange(AnimalID)
   
 nanimal <- nrow(trajectory)
 qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
 col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
 if (nanimal < length(col_vector)) colors <- sample(col_vector, nanimal) else colors <- sample(col_vector, nanimal,
                                                                                                     replace=TRUE)
 
 or.table <- sf.dat  %>%  as_tibble()  %>%  dplyr::select(-geometry)
 points.a <- left_join(points,or.table, join_by(AnimalID)) # join the original attributes back, doesn't completely solve
 
  #mapview(trajectory,zcol="AnimalID",color=colors)+mapview(points.a,zcol="AnimalID",cex=3,col.regions=colors,alpha.regions=1,legend=FALSE)
  mapview(trajectory,zcol="AnimalID",color=colors)+mapview(sf.dat,zcol="AnimalID",cex=3,col.regions=colors,alpha.regions=1,legend=FALSE)
  
}

# Simple and fast function Use MCP area to look for mortality like point patterns in the data #
#. inputs: projected sf object (not lat/lon), area threshold to filter results (m2)
#   output: mcp polygons with area column in m^2
FindMortPattern_byMCP <- function(sf.dat,area.threshold){
  
  # need at least 5 relocations to get MCP (caution with using dplyr group_by on sf objects)
  loc.per.animal <- sf.dat %>% st_drop_geometry() %>% group_by(AnimalID) %>%
    summarize(Fix_count =  n())
  drops <- which(loc.per.animal$Fix_count < 5,arr.ind=TRUE)
  excludes <- loc.per.animal$AnimalID[drops]
  if (length(drops)>0) sf.dat <- sf.dat[!sf.dat$AnimalID %in% excludes,]
  
  # convert sf to sp
  dat.sp <- sf.dat %>% select(AnimalID) %>% as("Spatial")
  
  mcps <- mcp(dat.sp,percent=100,unout="m2") # 95 allows an outlier
  
  mcps.sf <- st_as_sf(mcps)
  names(mcps.sf)[1] <- "AnimalID"
  MCP <- mcps.sf %>% filter(area < area.threshold)
  PointsOut <- sf.dat %>% filter(AnimalID %in% MCP$AnimalID)
  
  return(list(MCP=MCP,PointsOut=PointsOut))
  
}

# Alternative search for mort using dispersion msmts, see sfdep pkg
FindMortPattern_bySDD <- function(sf.dat,dist.threshold){
  
  # probably need at least 5 relocations  (caution with using dplyr group_by on sf objects)
  loc.per.animal <- sf.dat %>% st_drop_geometry() %>% group_by(AnimalID) %>%
    summarize(Fix_count =  n())
  drops <- which(loc.per.animal$Fix_count < 3,arr.ind=TRUE)
  excludes <- loc.per.animal$AnimalID[drops]
  if (length(drops)>0) sf.dat <- sf.dat[!sf.dat$AnimalID %in% excludes,]
  
  # compute std_distance measure to each point pattern by animalID
  dat.sdd <- sf.dat %>% group_by(AnimalID) %>% group_map(~ std_distance(.x))
  AnimalID <- sf.dat %>% select(AnimalID) %>% st_drop_geometry() %>% unique()
  sdd <- data.frame(AnimalID=AnimalID, SDD=unlist(dat.sdd))
  out.sdd <- sdd %>% filter(SDD < dist.threshold)
  
  return(out.sdd)
  
}

# A function to query a DB and return the most recent N GPS fixes for each animal unless a vector of IDs is provided
#   Inputs: tbl_db - a reference to the GPS table in the DB from dplyr function tbl
#           n   - number of records to return
#           ids - options vector of IDs to match in AnimalID column
FetchLastNFixes <- function(tbl_db, n, ids=NULL) {
  
  # query for appropriate data
  if (length(ids) > 0) {
  gps <- tbl_db %>% filter(AnimalID %in% ids) %>% group_by(AnimalID) %>% slice_max(order_by="acquisitiontime",n=n) %>% collect() 
  } else {
    gps <- tbl_db %>% group_by(AnimalID) %>% distinct() %>% slice_max(order_by="acquisitiontime",n=n) %>% collect()
  }
}

# a simple function to calculate prevalence in Movi data, ie. counting 'Detected' in data and dividing by total
calcPrevalence <- function(data) {
  
  n_detected <- length(which(data=="Detected", arr.ind=TRUE))
  prev <- n_detected/length(which(!(is.na(data)), arr.ind=TRUE))
  
  return(prev)
}
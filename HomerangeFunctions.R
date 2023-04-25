
# this script houses some functions to compute homeranges from GPS tables from our DB#

# required packages #
library(sf)
library(rgdal)
library(spatial)
library(adehabitatHR)
library(sp)
library(dplyr)
# required packages #


calculateHomerange <- function(gps, min.fixes, output.proj){
  
  # check for missing lat/lon and drop #
    missing.loc <- which(is.na(gps$latitude),arr.ind=TRUE)
    gps <- gps[-(missing.loc),]
  
  # check number of locations per animal, adehabitat requires >5/per
    loc.per.animal <- gps %>% group_by(AnimalID) %>%
      summarize(Fix_count =  n())
    drops <- which(loc.per.animal$Fix_count < min.fixes,arr.ind=TRUE)
    excludes <- loc.per.animal$AnimalID[drops]
    if (length(drops)>0) gps <- gps[!gps$AnimalID %in% excludes,]
  
    # ensure AnimalID is a factor
      gps$AnimalID <- factor(gps$AnimalID)
    
    # add coord ref for each list element crs= 4326 parameter assigns a WGS84 coordinate system #
     gps.sf <- st_as_sf(gps, coords = c("longitude", "latitude"), crs = 4326)
    
    # convert from sf to spatial points data frame (reqd by adehabitatHR)
      gps.sp <- gps.sf %>% select(AnimalID) # Keep just the "AnimalID" column of spatial object (again a restraint of the package)
      gps.sp <- as(gps.sf, "Spatial")
    
    
      # transform each data set in list from WGS84 to projected coords
      # and compute UD and polygon HR
        gps.collarIDT <- spTransform(gps.sp, output.proj)
        kud <- kernelUD(gps.collarIDT, grid=400,same4all = TRUE, h="href",extent=0.5)
        homeranges <- getverticeshr(kud)
        
      # construct output
        output <- list(homeranges=homeranges, UD=kud)
        
    return(output)
      
    
  
}
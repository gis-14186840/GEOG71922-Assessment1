# GEOG71922: Assessment1
# Author [14186840]

# 1.Setup and loading data

#set the working directory
setwd('E:/Manchester/S2GEOG71922_SE/Assessment1/Code/GEOG71922-Assessment1')

#load the required libraries
library(terra)
library(sf)  #terra and sf handling spatial data
library(dismo)  #dismo for downloading species data and evaluating models
library(glmnet)
library(maxnet)  #glmnet and maxnet for running Maxent models
library(precrec)  #precrec for model evaluation
library(dplyr)  #dplyr for data manipulation

#load in the land cover map
LCM=rast("LCMUK.tif")

#load in the Eurasian badger (Meles meles) data
meles=read.csv("Melesmeles.csv")

#Data cleaning
#subset the data to only include points with complete coordinates
meles<-meles[!is.na(meles$Latitude),]

#subset to cases with coordinate uncertainty less 1000 meters
meles<-meles[meles$Coordinate.uncertainty_m <= 1000,]

#Make spatial points layer
#create crs object
meles.latlong=data.frame(x=meles$Longitude, y=meles$Latitude)

#use coordinates object to create our spatial points object
meles.sp=st_as_sf(meles.latlong, coords=c("x", "y"), crs="epsg:4326")

# 2.Spatial cropping and projection

#set the extent to something workable
studyExtent<-ext(-4.2,-2.7,56.5,57.5) #list coordinates: min x, max x, min y, max y

#crop points to this area
meles.sp.crop<-st_crop(meles.sp, st_bbox(studyExtent))

#set the points to the same projection as the LCM layer
melesFin<-st_transform(meles.sp.crop,crs(LCM))

#crop the land cover data to the extent of the points data (plus 5km for buffers)
melesCoords<-st_coordinates(melesFin)
x.min <- min(melesCoords[,1]) - 5000
x.max <- max(melesCoords[,1]) + 5000
y.min <- min(melesCoords[,2]) - 5000
y.max <- max(melesCoords[,2]) + 5000
extent.new <- ext(x.min, x.max, y.min, y.max)

#crop the LCM raster to the extent
LCM_crop <- crop(LCM$LCMUK_1, extent.new)

#test drawing
plot(LCM)
plot(melesFin,add=TRUE)


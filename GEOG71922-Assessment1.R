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

#subset to cases with coordinate uncertainty <100 meters
meles<-meles[meles$Coordinate.uncertainty_m<100, ]

#Make spatial points layer
#create crs object
meles.latlong=data.frame(x=meles$Longitude, y=meles$Latitude)

#use coordinates object to create our spatial points object
meles.sp=st_as_sf(meles.latlong, coords=c("x", "y"), crs="epsg:4326")
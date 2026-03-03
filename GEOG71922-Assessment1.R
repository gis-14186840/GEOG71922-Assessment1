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

# 3.Prepare covariates

#access levels of the raster by treating them as categorical data
LCM_factor<-as.factor(LCM_crop)
lcm_classes<-levels(LCM_factor)[[1]][,1] 

#create an vector object called reclass
reclass_wood = ifelse(lcm_classes==1,1,0)

#combine with the LCM categories into a matrix of old and new values
RCmatrix=apply(cbind(lcm_classes, reclass_wood),2,as.numeric)

#asssign new values to LCM with reclassification matrix
broadleaf=classify(LCM_crop, RCmatrix)

#create an vector object called reclass Urban which is zero for all classes except tghe two urban classes in the LCM
urban=classify(LCM_crop, rbind(c(20,1), c(21,1)), others=0)

#aggregate LCM raster
broadleaf_agg=aggregate(broadleaf,fact=4,fun="modal")
urban_agg=aggregate(urban,fact=4,fun="modal")

#calculate the proportion of woodland within an 1800m circular neighborhood
wood_1800=focal(broadleaf_agg,w=focalMat(broadleaf_agg,1800,"circle"),na.rm=TRUE)

#calculate the proportion of urban areas within a 2300m circular neighborhood
urban_2300=focal(urban_agg,w=focalMat(urban_agg,2300,"circle"),na.rm=TRUE)

#stack the covariate layers together
allEnv=c(wood_1800, urban_2300)
names(allEnv)=c("broadleaf","urban")

#temp check result
#print(allEnv)

# 4.Generate background seeds and dataset

#create background points
set.seed(11)

#sample background - one point for every cell (9775)
back = spatSample(allEnv,size=2000,as.points=TRUE,method="random",na.rm=TRUE) 
back=back[!is.na(back$broadleaf),]
back=st_as_sf(back,crs="EPSG:27700")

#get environmental covariates at presence locations
eP=terra::extract(allEnv, vect(melesFin))

#bind together the presence data
Pres.cov=st_as_sf(cbind(eP,melesFin))
Pres.cov$Pres=1

#Remove the first column which is just an ID field.
Pres.cov=Pres.cov[,-1]

#get coordinates for spatial cross-validation later
coordsPres=st_coordinates(Pres.cov)

#drop geometry column
Back.cov=st_as_sf(data.frame(back,Pres=0))

#get coordinates of background points for cross validation later
coordsBack=st_coordinates(back)

#combine
coords=data.frame(rbind(coordsPres,coordsBack))

#assign coumn names
colnames(coords)=c("x","y")

#combine pres and background
all.cov=rbind(Pres.cov,Back.cov)

#add coordinates
all.cov=cbind(all.cov,coords)

#remove any NAs
all.cov=na.omit(all.cov)
all.cov=st_drop_geometry(all.cov)

#test result
print(all.cov)
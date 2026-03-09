# GEOG71922: Assessment1
# Author [14186840]

# 1.Setup and loading data

#set the working directory
setwd('E:/Manchester/S2GEOG71922_SE/Assessment1/Code/GEOG71922-Assessment1')

#load the required libraries
library(terra)
library(sf)  #terra and sf handling spatial data
library(dismo)  #dismo for downloading species data and evaluating models
library(maxnet)  #glmnet and maxnet for running Maxent models
library(dplyr)  #dplyr for data manipulation

#load in the land cover map
LCM=rast("LCMUK.tif")

#load in the Eurasian badger (Meles meles) data
meles=read.csv("Melesmeles.csv")

#data cleaning: filter complete coordinates, high precision, confirmed records
meles<-subset(meles,!is.na(meles$Latitude)&
                    Coordinate.uncertainty_m<=1000&
                    Identification.verification.status!="Unconfirmed")

#create spatial points object using WGS84
meles.sp=st_as_sf(meles,coords=c("Longitude","Latitude"),crs="epsg:4326")

# 2.Spatial cropping and projection

#set the extent to the study area
scot=st_read('scotSamp.shp')

#project squirrel data
melesFin.sp=st_transform(meles.sp, crs(LCM))

#crop points to the study area
melesFin=melesFin.sp[scot,]

#crop to the extent of the study area plus a little more
LCM_crop=crop(LCM$LCMUK_1,st_buffer(scot,dist= 1000))

#mask the LCM to this boundary
LCM_crop=mask(LCM_crop, scot)

# 3.Prepare covariates

#asssign new values to LCM with reclassification matrix
broadleaf=classify(LCM_crop,matrix(c(1,1),ncol=2),others=0)

#create an vector object called reclass Urban which is zero for all classes except tghe two urban classes in the LCM
urban=classify(LCM_crop,matrix(c(20,1,21,1),ncol=2,byrow=TRUE),others=0)

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

# 4.Generate background seeds and dataset

#create background points
set.seed(11)

#sample background points
back=spatSample(allEnv,size=2000,xy=TRUE,method="random",na.rm=TRUE) 
back$Pres=0

#extract coordinate and environmental covariate values
pres_coords=st_coordinates(melesFin)
eP=terra::extract(allEnv,pres_coords)

#combine the extracted data into dataframe
pres=data.frame(x=pres_coords[,1],y=pres_coords[,2], 
                  broadleaf=eP$broadleaf,urban=eP$urban,Pres=1)

#combine presence and background into dataframe
all.cov=na.omit(rbind(pres,back))

# 5.Model fitting: GLM vs Maxnet

#build the binomal glm model
glm_model=glm(Pres~poly(broadleaf,2)+poly(urban,2),family=binomial(link='logit'),
              data=all.cov)

#build the maxnet model
env_data=all.cov[,c("broadleaf","urban")]
maxnet_mod=maxnet(p=all.cov$Pres,data=env_data,classes="lq")

# 6.Test (evaluate) the model

#set number of folds to use
folds=5

#partiction presence and absence data according to folds using the kfold() function
kfold_pres<-kfold(all.cov[all.cov$Pres==1,],folds)
kfold_back<-kfold(all.cov[all.cov$Pres==0,],folds)

#create empty numeric vector to hold results
auc_glm<-numeric(folds)
auc_max<-numeric(folds)
opt_glm<-numeric(folds) 
opt_max<-numeric(folds)

#for loop to iterate over folds
for (i in 1:folds) {
  #for presence values, select all folds which are not 'i' to train the model
  #bind presence and background training data together
  datatrain<-rbind(all.cov[all.cov$Pres==1,][kfold_pres!=i, ], 
                   all.cov[all.cov$Pres==0,][kfold_back!=i,])
  
  #the remaining fold is used to test the model
  test_pres<-all.cov[all.cov$Pres==1,][kfold_pres==i,]
  test_back<-all.cov[all.cov$Pres==0,][kfold_back==i,]
  
  #glm model trained on presence and absence points
  glm_cv<-glm(Pres~poly(broadleaf,2)+poly(urban,2),family=binomial,data=datatrain)
  
  #use testing data for model evaluation
  pred_glm.p<-predict(glm_cv,test_pres,type= "response")
  pred_glm.a<-predict(glm_cv,test_back,type="response")
  eval_glm<-evaluate(p=pred_glm.p,a=pred_glm.a)
  auc_glm[i]<-eval_glm@auc
  
  #list and select corresponding values for TPR(true positive rate)and TNR(true negative rate)are highest
  opt_glm[i]<-eval_glm@t[which.max(eval_glm@TPR +eval_glm@TNR)]
  
  #maxnet evaluation
  max_cv<-maxnet(p=datatrain$Pres,data=datatrain[,c("broadleaf","urban")])

  #use dismo::evaluate to predict AUC and threshold
  pred_max.p<-as.numeric(predict(max_cv, test_pres[,c("broadleaf", "urban")],type ="cloglog"))
  pred_max.a<-as.numeric(predict(max_cv, test_back[,c("broadleaf", "urban")],type ="cloglog"))
  eval_max_dismo<-evaluate(p=pred_max.p, a=pred_max.a)
  
  #extract the AUC and max threshold
  auc_max[i]<-eval_max_dismo@auc
  opt_max[i]<-eval_max_dismo@t[which.max(eval_max_dismo@TPR+ eval_max_dismo@TNR)] }

#print results
print(paste("GLM 5-Fold Mean AUC:",mean(auc_glm)))
print(paste("Maxnet 5-Fold Mean AUC:",mean(auc_max)))
print(paste("GLM Optimal Threshold:",mean(opt_glm)))
print(paste("Maxnet Optimal Threshold:",mean(opt_max)))

# 7.Prediction and mapping

#set plotting scheme
par(mfrow=c(1, 1))

#Plot the study area and basemap (1800m broadleaf woodland proportion)
plot(wood_1800, main="Study Area & Occurrence Data", 
     col=terrain.colors(100), axes=TRUE)

#add the 2000 generated background points. "pch" sets the symbol type and "cex" sites the size
plot(st_geometry(back_sf), col="grey", pch=20, cex=0.3, add=TRUE)

#Do the same for Eurasian badger presence records
plot(st_geometry(melesFin), col="red", pch=4, cex=0.8, lwd=1.5, add=TRUE)

#add legend
legend("bottomright", legend=c("Presence (M. meles)", "Background (Pseudo-absences).)"), 
       col=c("red", "grey"), pch=c(4, 20), pt.cex=c(0.8, 0.5), bg="white", cex=0.7)

#generate model response curves
#restore 1row, 2column plot layout
par(mfrow=c(1, 2))

#plot the response curve for Broadleaf woodland
wood_seq=seq(min(all.cov$broadleaf),max(all.cov$broadleaf),length.out=100)
test_data_wood=data.frame(broadleaf=wood_seq,urban=median(all.cov$urban))
pred_wood=predict(glm_model,newdata=test_data_wood,type="response")
plot(wood_seq,pred_wood,type="l",col="green",lwd=2,
     xlab="Broadleaf Woodland (1800m proportion)",ylab="Probability of Occurrence",
     main="Response to Woodland")

#plot the response curve for Urban density
urban_seq=seq(min(all.cov$urban),max(all.cov$urban),length.out=100)
test_data_urban=data.frame(broadleaf=median(all.cov$broadleaf),urban=urban_seq)
pred_urban=predict(glm_model, newdata=test_data_urban,type="response")
plot(urban_seq,pred_urban,type="l",col="red",lwd=2,
     xlab="Urban Density (2300m proportion)",ylab="Probability of Occurrence",
     main="Response to Urbanization")

#predict and inspect the output
map_glm=predict(allEnv,glm_model,type="response",na.rm=TRUE)
map_maxnet=predict(allEnv, maxnet_mod,type="cloglog",clamp=FALSE,na.rm=TRUE)

#make comparison plot
par(mfrow=c(2, 2))
plot(map_glm,main="GLM Continuous Probability")
plot(map_maxnet,main="Maxent Continuous Probability")
plot(map_glm>mean(opt_glm),main=paste("GLM Binary(T>", round(thresh_glm,3),")"),col=c("lightgrey","darkgreen"))
plot(map_maxnet>mean(opt_max),main=paste("Maxent Binary(T>",round(thresh_max,3),")"),col=c("lightgrey","darkgreen"))

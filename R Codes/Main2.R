# set the working directory
w_dir <- "/home/hbg/Desktop/Spain_american_mink_project"
setwd(w_dir)

# imported important libraries to handle geo-spatio-temporal data
library(spatstat)
library(raster)
library(terra)
library(sp)

map_folder<- paste(w_dir,"/maps/",sep="")

# Read occurrence only data of American mink population in Spain
occ = read.csv(paste(w_dir,"/data/Mink_spain_occ.csv",sep=""))
str(occ)


hist(occ$year)  # got estimates of no. of data points each year

occ <- occ[!(occ$year == 2020),]                        # removed 2020 as it has very less recorded data points

# removed all the points without values in geo-spatial or year column
occ<- occ[!is.na(occ$year),]
occ<- occ[!is.na(occ$longitude),]
occ<- occ[!is.na(occ$latitude),]


# converted Occurrence points to spatial points dataframe with CRS WGS84 datum
occ.points <- SpatialPoints(occ[,c("longitude", "latitude")], proj4string = CRS("EPSG:4326"))    
number_of_points  <- length(occ.points)
n <- number_of_points
occ.datapts = SpatialPointsDataFrame(occ.points, data.frame(ID=1:n,state=1,year=as.numeric(occ$year)))

Spain <- shapefile(paste(map_folder,"Spain_vector/","Spain.shp",sep=""))       # Shape file of Spain national boundary for the reference of visualizations



occ.datapts <- occ.datapts[Spain,]                                     # excluded all the outlier points mistakenly reported outside Spain boundary
str(occ.datapts)

op <- par(mfrow=c(1,1))                                                 #plot_settings (1*1) figures



# First map in folder Res_tiffs is predicted occurrence probability of American mink (0 to 1)  in Spain environment using maxent algorithm 
# MaxEnt contrasts the conditions between presence locations and the study area to estimate a presence probability surface. 
# central idea in MaxEnt is to search for a probability distribution having a maximum entropy.
# Used Maxent software for modeling species niches and distributions (Version 3.4.1) [Steven J. Phillips, Miroslav Dudík, Robert E. Schapire]

# other two maps in folder are restricted area where mink less likely to be found on the basis ecological restriction of elevation > 1500 m and river distance > 60 Km 
# these two maps generated DEM map and River shape_file with the help of Raster calculator method of Q-GIS

res.tiffs <- list.files(paste(map_folder,"Res_tiffs",sep=""), full.names=TRUE)          
basename(res.tiffs)

r1<- raster(res.tiffs[1])
max.ml.p2 <- as(r1, "SpatialGridDataFrame")

# Generated 10% psuedo-absence points from area with less than 1% probability of finding American mink

max.ml.p2$absence = ifelse((max.ml.p2[[1]]<0.01), 1 , NA )                           
dens.var <- spatstat.geom::as.im(sp::as.image.SpatialGridDataFrame(max.ml.p2["absence"]))
plot(r1)
pnts.new1 <- rpoint((number_of_points * 0.10), f=dens.var)
points(pnts.new1, col="red")



r2 <- stack(res.tiffs[2],res.tiffs[3])
max.ml.p2 <- as(r2, "SpatialGridDataFrame")

# Random 5% Psuedo-absence points from available restricted area
max.ml.p2$absence = ifelse((max.ml.p2[[1]]==1) | (max.ml.p2[[2]]==1), 1 , NA ) 
dens.var <- spatstat.geom::as.im(sp::as.image.SpatialGridDataFrame(max.ml.p2["absence"]))
plot(r2[[1]]|r2[[2]])
pnts.new2 <- rpoint((number_of_points * 0.05), f=dens.var)
points(pnts.new2, col="red")

# Combined 15% Psuedo-absence points 
pnts.new <- superimpose(pnts.new1,pnts.new2)
plot(r1)
points(pnts.new, col="red")


pnts.new <- as(pnts.new, "SpatialPoints")
proj4string(pnts.new) <- CRS("+init=epsg:4326")

# Created spatial data-frame of all the psuedo-absence points by assigning random year with probability of that year being present in Occurrence dataset
n <- number_of_points
abs.datapts = SpatialPointsDataFrame(pnts.new, data.frame(ID=(n+1) : (n+length(pnts.new)), state=0, year= sample(occ.datapts$year,length(pnts.new))))
summary(abs.datapts)


#library(maptools)

occ_abs.matrix <-  occ.datapts + abs.datapts # Combined occurrence and absence data set

f2 <- paste(w_dir,"/data/",sep="")

write.csv(occ_abs.matrix,paste(f2, "occ_abs_A-mink.csv", sep=""))       # Saved Combined occurrence and absence data set for future use
str(occ_abs.matrix)
summary(occ_abs.matrix)

# initialized empty lists to store whole Spain environmental covariates data-sets and mink datasets binded with these covariates
g4km_l = c()
occ.final_l = c()


year_list <- unique(occ_abs.matrix$year)
print(year_list)                    # List of all unique years in Occurrence dataset

# Loop to load Spain environmental covariates data from geo-tiff mpas (5Km Resolution) and bind them to mink occurence_absence datasets
for (i in year_list){
year <- i
occ.new <- occ_abs.matrix[occ_abs.matrix$year == year ,]
occ.new$year= year

# Loading all the Dynamic covariates that changes every year
# List are - skin_temp of earth surface, surface solar radiation, Snow coverage, Corine_land_coverage, chelsa bio and leaf area index in high and low vegetation
# Each geo-tiffs have two bands one for jan (coldest month in Spain) and one for July (hottest month in Spain)

print(paste(w_dir,"/maps/Dynamic/",as.character(year),sep =""))
eco.tifs1 = list.files(paste(w_dir,"/maps/Dynamic/",as.character(year),sep=""), glob2rx("*.tif$"), full.names=TRUE)
basename(eco.tifs1)

# Loading all the Static covariates that remains constant every year
# List are - Elevation, Distance_from_water and Water_wetness
eco.tifs2 = list.files(paste(w_dir,"/maps/Static/",sep=""), glob2rx("*.tif$"), full.names=TRUE)
basename(eco.tifs2)

eco.tifs <- c(eco.tifs1 , eco.tifs2)
basename(eco.tifs)


# stacking all the covariates maps of single year into one stack
st <- raster::stack(eco.tifs)

# converted the stacked raster into pixel-dataframe  of covariates
g4km = st
g4km = as(g4km, "SpatialGridDataFrame")
g4km = as(g4km, "SpatialPixelsDataFrame")

g4km_l = append(g4km_l,g4km)  # appended each year the covariates dataframe into a list to further use for predictions

# find the associated value of covariates of occurrence/absence points by overlapping and binding it with covariates data frame
occ.regmat <- sp::over(occ.new, g4km)
occ.final <- cbind(occ.new, occ.regmat)

occ.final_l = append(occ.final_l,occ.final)   # appended combined each year dataframe into a list to further use for training
}

# joining every year dataset into one for training the model
occ.final <- occ.final_l[[1]]

for (i in c(2:length(year_list))){
occ.final <-  occ.final + occ.final_l[[i]] 
}
occ.final$ID <- (1:length(occ.final))

summary(occ.final)



# Splitted the original dataset into 80-20 Train-Test dataset
require(caTools)
set.seed(151) #101 #151
sample = sample.split(occ.final$state, SplitRatio = .80)
train = subset(occ.final, sample == TRUE)
test  = subset(occ.final, sample == FALSE)

str(train)

# train_new <- train[train$year == 2010, ]
# train_new <- train_new[c(sample(1:length(train_new$year), 400)),]
# 
# year_list <- unique(train$year)
# 
# for (i in year_list[-1]){
#   train_y <- train[train$year == i, ]
#   train_new <- train_new + train_y[c(sample(1:length(train_y$year), 222)),]
# }

# visualizing the final training dateset

plot(Spain)
points(train[train$state == 1,], col ="blue")
points(train[train$state == 0,], col ="red")

########################################################################################
######## EMl training and prediction ########################

# imported all the important libraries required for Ensemble-machine Learning (EML)
library(ranger)
library(rgdal)
library(geoR)

library(glmnet)
library(xgboost)
library(kernlab)

library(deepnet)
library(forestError)
library(mlr)
library(dplyr)
library(stats)


# Formula created for the machine learning by encompassing all the covariates that can affect state of occurrence
covv = paste(names(train[-c(1:3)]), collapse="+")
str(covv)
fm.fs = stats::as.formula(paste("state ~ ", paste("year+",covv)))
str(all.vars(fm.fs))


# Imported other R file created with all the functions related to Ensemble Machine Learning
source("R Codes/Eml.R")

#converted data_type of state into factor that can be used for classification Learning
train$state <- as.factor(train$state)
train$year <- as.numeric(train$year)
str(train)

f2 <- paste(w_dir,"/Output/",sep="")

# tuned the parameters using training data in according to use it for  Ensemble of classif.ranger, classif.xgboost, classif.glmnet learners and removed unwanted covariates
# ranger is based on random forest and glmnet is baded on neuaral-net both model had shown importance in geo-spatial datasets
# xgboost is gradient-boosting decision trees is good for fine-tune.
# All these function is written to use parallel computation of all the CPU cores.
tnd.ml = tune_learners(data = train@data , formula = fm.fs , 
                       blocking = factor(train$ID), out.dir=f2)

# Run the tuned super_learner of Ensemble of classif.ranger, classif.xgboost, classif.glmnet learners
t.m = train_sp_eml(data = train@data, tune_result = tnd.ml, 
                   blocking = as.factor(train$ID), out.dir = f2)

summary(t.m$learner.model$super.model$learner.model)

# in the Summary of EML model Estimate is weightage of different models and stars is significance

# To get the 10 most significant covariates factor for the training
xl <- as.data.frame(mlr::getFeatureImportance(t.m[["learner.model"]][["base.models"]][[1]])$res)
xl$relative_importance = round(100*xl$importance/sum(xl$importance), 1)
xl = xl[order(xl$relative_importance, decreasing = T),]
xl$variable = paste0(c(1:length(t.m$features)), ". ", xl$variable)

xl[1:10,]


# imported the Grid of 5Km of whole Spain Country for prediction boundary.
Spain.grid <- raster::shapefile("maps/Spain.grid/Spain_grid_0,05.shp")

Spain.grid_005 <- SpatialPixelsDataFrame(Spain.grid, data.frame(ID=1:(length(Spain.grid))))


# Spain.grid <- raster::shapefile("maps/Spain.grid/Spain_grid_0,01.shp")
# 
# str(Spain.grid)

#Spain.grid_001 <- SpatialPixelsDataFrame(Spain.grid, data.frame(ID=1:(length(Spain.grid))))



# Loop for predicting the occurrence probability using EML fitted model each year 
for (i in c(1:length(year_list))){

# Dataframe of covariates over grid of 5Km to predict the occurrence probability  
pred.cov <- sp::over(Spain.grid_005, g4km_l[[i]])
pred.cov <-  cbind(Spain.grid_005, pred.cov)

str(pred.cov)

multiplier=100
min.md=1

out.year <- year_list[i]
pred.cov$year <- out.year

# Removed all the Rows having any NA values of Covariates
for (j in c(1:length(pred.cov@data))){
  pred.cov<- pred.cov[!is.na(pred.cov@data[[j]]),]
}

cc = complete.cases(pred.cov@data)
g1km = pred.cov
g1kx = pred.cov[1]


library("RColorBrewer") # for better colour-scheme of plots

#predict states using trained Super learner and grid-covariates dataframe
pred <- predict(t.m, newdata=g1km@data[,t.m$features])$data$prob.1 
g1kx@data[cc,"pred"] = round(pred * 100)
str(g1kx)
# Saved the predicted map into geo-tiff file 
rgdal::writeGDAL(g1kx["pred"], paste0(f2, out.year, "_f.tif"), type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")

#visualizing the predicted plot
op <- par(mfrow=c(1,2))
plot(raster(g1kx["pred"]), col=brewer.pal(n = 9, name = "GnBu"),
     main="Prediction EML", axes=FALSE, box=FALSE)

#points(occ.datapts[occ.datapts$year == out.year, ], pch="+", cex=.8)

# extracted the allowed error for prediction from learned model
out.p <- as.matrix(as.data.frame(mlr::getStackedBaseLearnerPredictions(t.m, newdata=g1km@data[,t.m$features])))
g1kx@data[cc,"model.error"] = matrixStats::rowSds(out.p*100, na.rm=TRUE)
g1kx$model.error <- ifelse(g1kx$model.error<min.md, min.md, g1kx$model.error)

# Saved the Permissible prediction error map into geo-tiff file 
rgdal::writeGDAL(g1kx["model.error"],paste0(f2, out.year, "error.tif"), type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")

#visualizing the prediction error plot
plot(raster(g1kx["model.error"]), col=rev(bpy.colors()),
     main="Permissible Prediction errors", axes=FALSE, box=FALSE)
#points(occ.datapts[occ.datapts$year == out.year, ], pch="+", cex=.8)
}


par( mfrow= c(1,1) , mar= c(1,2,5,1)) 


# use test dataset to calculate model relative accuracy
str(test)

cc = complete.cases(test@data)
pred_test_m = test[1]

for (j in c(1:length(test@data))){
  test<- test[!is.na(test@data[[j]]),]
}

pred_test <- predict(t.m, newdata=test@data[,t.m$features])$data$prob.1 
pred_test_m@data[cc,"pred"] = pred_test

str(pred_test_m)
pred_test_m$error <- 1


pred_test_m$error <- abs(pred_test_m$pred-test$state)

pred_test_m<- pred_test_m[!is.na(pred_test_m$error),]

# how much prediction probability deviated from original value 1 and 0
Avg_error <- sum(pred_test_m$error)*100/length(pred_test_m$error)
Avg_error

pred_test_m$state <- 1
for (i in c(1:length(pred_test_m$state))){
  if(pred_test_m$pred[i] < 0.5){
    pred_test_m$state[i] <- 0
  }
}

conv_mat <- table( pred_test_m$state == test$state) # create a table of Matching of states of test and predicted

Accuracy <- conv_mat[[2]]/ ( conv_mat[[2]] + conv_mat[[1]])    # True Matching/ ( True + False)
Accuracy


library(plotKML)
# # Load and animate all the prediction maps from 2010-2019

es.tifs = list.files(f2, glob2rx("*_f.tif"), full.names = TRUE)
basename(es.tifs)
spain1km = raster::brick(raster::stack(es.tifs))

# animate(spain1km, pause = 0.25, n=1)


# Library to generate trend map using temporal trend of each pixel from 2010-2019
library(greenbrown)

data(SAGA_pal)
trendmap <- TrendRaster(spain1km, start=c(2010, 1), freq=1, breaks=1) 

par( mfrow= c(1,1) , mar= c(1,2,5,1))

plot(trendmap[["SlopeSEG1"]], 
     col=rev(SAGA_pal[["SG_COLORS_GREEN_GREY_RED"]]), 
     zlim=c(-1.5,1.5), main="Slope SEG1")



# Generated another trend Map manually using temporal trend of each pixel of predicted raster for every year
xs = as(spain1km, "SpatialGridDataFrame")
in.years = year_list

str(in.years)
# used parallel computing using multiple cores of CPU
cl <- parallel::makeCluster(parallel::detectCores())
parallel::clusterExport(cl, c("in.years"))
# select ONLY pixels that change in time
sd.pix = which(!parallel::parApply(cl, xs@data, 1, sd, na.rm=TRUE)==0)
# derive beta per pixel:
# Because probablities are binary variable used logit transformation
betas = unlist(parallel::parApply(cl, xs@data[sd.pix,], 1, function(i) { 
  try( round( coef(lm(y~x, data.frame(x=in.years, 
                                      y=boot::logit(ifelse(as.vector(i)<=0, 0.1, ifelse(as.vector(i)>=100, 99.9, as.vector(i)))/100))))[2] *1000 ) 
  ) }))
# Saved it to geo-tiff raster file
spain.trend = xs[1]
spain.trend@data[,1] = 0
spain.trend@data[sd.pix,1] = as.numeric(betas)
rgdal::writeGDAL(spain.trend[1], paste(f2,"spain_trend_aedes_1km.tif"), type="Int16", 
                 mvFlag=-32768, options=c("COMPRESS=DEFLATE"))
parallel::stopCluster(cl)

plot(spain.trend[1], 
     col=rev(SAGA_pal[["SG_COLORS_GREEN_GREY_RED"]]), main="Slope SEG1")

print(paste("Avergae error : ", as.character(Avg_error)))
print(paste("Avergae error : ", as.character(Accuracy)))

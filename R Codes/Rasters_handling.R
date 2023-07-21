# this is the code to convert Raster maps into ascii formatto use them in Maxent software and Rangeshifter Algorithm
# set the working directory
w_dir <- "/home/hbg/Desktop/Spain_american_mink_project"
setwd(w_dir)

f<- paste(w_dir,"/maxENT/ENV",sep="")

eco.tifs <- list.files(paste(f), glob2rx("*.tif$"), full.names=TRUE)
basename(eco.tifs)

library(raster)

r <- stack(eco.tifs)

plot(r[[-5]])

for(i in 1:nlayers(r)){
  
band<-r[[i]] 
writeRaster(band,paste(f,"asc2/",as.character(names(r[[i]])),'.asc', sep='', format="ascii")) 
}



f<- "/home/hbg/Desktop/Mink spain (copy)/imp maps/Range_shift/"

eco.tifs <- list.files(paste(f), glob2rx("*.tif$"), full.names=TRUE)
basename(eco.tifs)

library(raster)
r<-c()
r[[1]] <- stack(eco.tifs[[1]])
r[[2]] <- stack(eco.tifs[[2]])
#plot(r)

for(i in 1:2){
  
  band<-r[[i]] 
  writeRaster(band,paste(f,"asc/",as.character(names(r[[i]])),'.asc', sep='', format="ascii")) 
}

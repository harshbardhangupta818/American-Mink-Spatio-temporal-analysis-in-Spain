# set the working directory
w_dir <- "/home/hbg/Desktop/Spain_american_mink_project"
setwd(w_dir)

# Read occurrence only data of American mink population in Spain
occ = read.csv(paste(w_dir,"/data/Mink_spain_occ.csv",sep=""))
str(occ)

dates <- as.Date( occ$Date, "%m/%d/%Y")   #correct the format of date
occ$Date <- dates 


# removed all the points without values in geo-spatial or recorded date column
occ<- occ[!is.na(occ$Date),]
occ<- occ[!is.na(occ$longitude),]
occ<- occ[!is.na(occ$latitude),]

# imported important libraries to handle geo-spatio-temporal data
library(spacetime)
library(plotKML)
library(spatstat)
library(sp)

# converted Occurrence points to spatial-temporal dataframe with CRS WGS84 datum
sp_ST <- STIDF(SpatialPoints(occ[,c("longitude", "latitude")], proj4string = CRS("EPSG:4326")), 
               occ$Date, data.frame(individualCount=occ$data.Type))
data(SAGA_pal)


## Directly plot in Google Earth:
#plotKML(sp_ST, folder.name ="/home/hbg/Desktop/Mink spain (copy)/" , file.name = "American_Mink.kml" , colour_scale=SAGA_pal[[1]])

library(raster)
# Import Water Shape_file
water <- raster::shapefile("maps/KML/Water_ways.shp")

# open a KML file
kml_open("Plots/Mink_Label.kml")
# Add Spatiotemporal occurrence data layer
kml_layer(sp_ST, colour=sp_ST@data[,1], shape="http://maps.google.com/mapfiles/kml/paddle/blu-blank.png", size = 0.7,
          colour_scale=SAGA_pal[[1]], points_names="" ,  driver="LIBKML" , layer="multipoint")
# Add Water Shape-file layer
kml_layer(water, colour = "white",  driver="LIBKML", layer = "linestring")
kml_close("Plots/Mink_Label.kml")


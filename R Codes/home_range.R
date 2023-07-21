# set the working directory
w_dir <- "/home/hbg/Desktop/Spain_american_mink_project"
setwd(w_dir)

# Read occurrence only data of American mink population in Spain
occ = read.csv(paste(w_dir,"/data/Mink_spain_occ.csv",sep=""))
str(occ)


# imported important libraries to handle geo-spatio-temporal data
library(terra)
library(sp)
library(sf)

# library foir home range analysis
library(adehabitatHR)
library(dplyr)

# removed all the points without values in geo-spatial or year column
occ<- occ[!is.na(occ$year),]
occ<- occ[!is.na(occ$longitude),]
occ<- occ[!is.na(occ$latitude),]



#check if any year have less than 5 ( needed for convex polygon) data points using dplyr
check<- occ %>% group_by(year) %>% summarise(loop=length(year)) %>% dplyr::filter(loop<5)

#convert this to a dataframe
check<- as.data.frame(check)

#subset occ to exclude those individuals found to have less than 5 data points

occ<-occ %>% anti_join(check)


# convert Occurrence points to spatial points dataframe with CRS WGS84 datum
occ.points <- SpatialPoints(occ[,c("longitude", "latitude")], proj4string = CRS("EPSG:4326"))

# Transform the CRS (coordinate reference system) to local metric coordinate of Spain mainland with easting and northing 
occ.points <- spTransform(occ.points, CRS("EPSG:2062"))
n  <- length(occ.points)
occ.datapts = SpatialPointsDataFrame(occ.points, data.frame(ID=1:n,year=occ$year))                          

library(rgdal)

summary(occ.datapts)
# imported the Grid of 1Km of whole Spain Country for saving polygons.

Spain.grid <- raster::shapefile("maps/Spain.grid/Spain_grid_1km_UTM.shp")

str(Spain.grid)
Spain.points <- spTransform(Spain.grid, CRS("EPSG:2062"))
Spain.grid_1km <- SpatialPixelsDataFrame(Spain.points, data.frame(ID=1:(length(Spain.grid))))


# Parallel library for parallel computing to save time
library(parallel)

year_list <- unique(occ.datapts$year)                     # List of all unique years in Occurrence dataset
print(as.character(year_list[9:11]))

# initialized empty lists to store area of home range of corresponding year
year_range <-c()
home_range <-c()

numCores <- detectCores()
numCores

library(raster)                             #For rasters handelling

y_list = seq.int(1,length(year_list) , by = numCores)
#for (i in y_list){
# Feeded 4 year at time as num_cores is 4
#years <- as.character(year_list[i:i+3])
years <- as.character(year_list[9:11])
area_fx <- function(years) {
  year <- years
  year_range <- append(year_range,as.numeric(year))
  occ.new <- occ.datapts[occ.datapts$year <= year ,]     # home range will be polygon made from all the points in that year along with points the previous year
  occ.new$state =1
  str(occ.new@proj4string)
  str(Spain.grid_1km@proj4string)


  # Creating a local convex polygon with points not more than total distance of 1200 Km
  alc <- LoCoH.a(occ.new[,3], a=1200000, unin = "m",
                 unout =  "km2", duplicates="remove")
  
  # convert and save that polygon as raster 
  par(mfrow= c(1,1))
  uur <- MCHu.rast(alc, Spain.grid_1km, percent=95)
  rast <- raster(uur,1)
  plot(rast,)
  points(occ.new)
  f<- "Home_range/"
  writeRaster(rast, paste(f,year,"homerange.shp"), overwrite=TRUE)
  
  rm(alc)
  rm(rast)
  rm(uur)
  
  # calculate the area by creating a local convex polygon with points not more than total distance of 1200 Km encompassing 95% points leaving 5 % outliers
  ar <- LoCoH.a.area(occ.new[,3], arange = c(1200000) , percent = 95,unin = "m",
                     unout =  "km2", duplicates="remove")

  rs <- ar[[1]]
  rm(ar)
  result<- c(year,rs)
}

system.time({
results2 <- mclapply(years, area_fx, mc.cores = numCores)})
#}

# repeat the above process for all the years from year list 


year_range <- append(year_range, unlist(results2)[ c(TRUE,FALSE) ])
home_range <- append(home_range, unlist(results2)[ c(FALSE,TRUE) ])


print(year_range)
print(home_range)

df <- data.frame(year_range, home_range)

# save home range area trend
write.csv(df,paste("Home_range/home_range.csv"))

################################
# plot home range expansion trend
hr <- read.csv("Home_range/home_range.csv")
str(hr)
Spain_area <- 505990
perc_home_range <- hr$home_range/Spain_area *100
plot(hr$year_range, perc_home_range , xlab ="year", ylab ="percent area(%)")
lines(hr$year_range, perc_home_range)


# Fitted a sigmoidal curve on the expansion trend
library(minpack.lm)
df1 <- data.frame(x=hr$year_range, y =  perc_home_range )
G <- expression(L / (1 + exp(-(k)*(x-x0))) + b)


Mm <- function(x, L, x0,k, b) {eval(G)}

model= nlsLM(y~Mm(x,L,x0,k,b), data = df1, start= list(L= max(df1$y),  x0 = median(df1$x), k = 1, b=min(df1$y)), algorithm = "LM", control = list(maxiter = 1000))

s_ <-  seq(2010, 2020, 0.01)

funx <- function(model_p,x) {
  L <- summary(model_p)$parameters[1]
  x0 <- summary(model_p)$parameters[2]
  k <- summary(model_p)$parameters[3]
  b <- summary(model_p)$parameters[4]
  y = Mm(x,L,x0,k,b)
  return (y)
}


pred1= funx(model,s_)



plot(hr$year_range, perc_home_range , xlab ="year", ylab ="area(%)", title("Home Range expansion"))
lines(s_, pred1, lty = 3 , lwd = 2, col ="red")


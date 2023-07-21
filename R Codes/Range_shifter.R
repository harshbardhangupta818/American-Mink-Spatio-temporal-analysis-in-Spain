library(RangeShiftR)
library(raster)
library(RColorBrewer)
library(rasterVis)
library(latticeExtra)
library(viridis)
library(grid)
library(gridExtra)

# set the path for working directory:

dirpath = "/home/hbg/Desktop/Spain_american_mink_project/RangeShift/"
setwd(dirpath) 

dir.create("Inputs", showWarnings = TRUE)
dir.create("Outputs", showWarnings = TRUE)
dir.create("Output_Maps", showWarnings = TRUE)


# Loaded Mean EML probability distribution map and actual distribution map           
Spain.map <- raster("Inputs/EML_result_avg_prob_madrid.tif")
SpDist <- raster("Inputs/Occ_2019_5km.tif")

values(SpDist)[values(SpDist) != 1] <- NA      # set all other value than 1 of actual distribution  to Na for better visualisation

#### can skip these steps if maps already generated

summary(values(Spain.map))
m = 0.99                     # set max value from above summary


# convert the Env suitable probability map to 3 habitat category
values(Spain.map)[(values(Spain.map) >= 2/3*m)] <- 3
values(Spain.map)[(values(Spain.map) >= 1/3 *m & values(Spain.map) <2/3 *m)] <- 2
values(Spain.map)[(values(Spain.map) >= 0 & values(Spain.map) <1/3 *m)] <- 1


summary(values(Spain.map))

# Plot actual distribution over habitat Category map with
plot(Spain.map, col=brewer.pal(n = 3, name = "Spectral"), axes=F)
plot(rasterToPolygons(SpDist, dissolve=F), add=T)

# convert the habitat map into ascii format with 5km resolution
r<- Spain.map
r@file@nodatavalue <- -9999
res(r) <- c(5000,5000)
res(r)
str(r)

writeRaster(r, format="ascii", filename = paste0(dirpath, "Inputs/Spain_EML_map_mink_3cl"), datatype = 'INT1S', overwrite = T)

# convert the actual distribution map into ascii format with 5km resolution

r2<- SpDist
values(r2)[is.na(values(SpDist))] <- 0
r2@file@nodatavalue <- -9999
res(r2) <- c(5000,5000)
str(r2)
writeRaster(r2, format="ascii", datatype = 'INT1S', filename = paste0(dirpath,"Inputs/occ_dist_2019_5km"), overwrite = T)

length(values(r2)[values(r2) == 1])



###################################################################

carrycap <- c(0, 0.0012, 0.002)                 # carrying capacity of each habitat type indv per hectare ( 1 ind in 25 Km2 cell = 0.0004 in 1 ha)

# import landscape File and actual distribution file
land <- ImportedLandscape(LandscapeFile = "Spain_EML_map_mink_3cl.asc", 
                          Resolution = 5000, 
                          Nhabitats = 3, 
                          K_or_DensDep = carrycap, 
                          SpDistFile = "occ_dist_2010_5km.asc", 
                          SpDistResolution = 5000)

(trans_mat <- matrix(c(0, 0.9, 0.4, 0.75), nrow = 2, byrow = F))  #  Two stage juvenile stage, the adult stage
# fecundity- 0.6, with a juvenile survival probability - 0.9 in their first year and  Adult survival probability - 0.75

stg <- StageStructure(Stages=2,           # 1 juvenile + adult stage
                      TransMatrix=trans_mat, 
                      MaxAge=10, 
                      SurvSched=2, 
                      FecDensDep=T)

demo <- Demography(StageStruct = stg,
                   ReproductionType = 0) #  Asexual female model

# check it for Rmax removing stagestrc for 0.6, 0.8, 1.0, 1.2, 1.4

disp <-  Dispersal(Emigration = Emigration(DensDep=T, StageDep=T, 
                                           EmigProb = cbind(0:1,c(1.0,0),c(10.0,0),c(1.0,0)) ),  # only juvenile transfer with 100% probability
                   Transfer = DispersalKernel(Distances = 40000) ,  # unit is m, so Mean = 40Km
                   Settlement = Settlement() )


init <- Initialise(InitType = 1, # = initialisation from a loaded species distribution map
                   SpType = 0,   # = all suitable cells within all distribution presence cells
                   InitDens = 2,
                   IndsHaCell = 2,  # 2 ind per initial cell of 25km2
                   PropStages = c(0.0,1.0),
                   InitAge = 0) 

sim_0 <- Simulation(Simulation = 0, 
                    Replicates = 5, 
                    Years = 10,
                    OutIntPop = 1,                             #output each year
                    OutIntOcc = 1,
                    OutIntRange = 1)

s <- RSsim(land = land, demog = demo, dispersal = disp, simul = sim_0, init = init)

s


validateRSparams(s)

# Run the final simulation
RunRS(s, dirpath)$Errors

# read 'range' output into a data frame

range_df <- readRange(s, dirpath)

# plot trajectories of all individual runs and overlay with mean:
par(mfrow=c(1,2))
plotAbundance(range_df)
plotOccupancy(range_df)

par(mfrow=c(1,2))
plotAbundance(range_df, rep=F, sd=T)
plotOccupancy(range_df, rep=F, sd=T)


# read pop output table into a data frame

pop_df <- read.table("Outputs/Batch1_Sim0_Land1_Pop.txt", header = TRUE)
extent(Spain.map)@xmin
str(pop_df)

# convert relative x,y coordinates of 5km resolution to original coordinates with CRS (ESPG:2062)
pop_df$x <- (pop_df$x * 5000) + extent(Spain.map)@xmin
pop_df$y <- (pop_df$y * 5000) + extent(Spain.map)@ymin
str(pop_df)


stack_pop <- function(pop_df, ext, rep=NULL, mask=NULL){
  # This function takes the population data frame output from RangeShiftR and turns it into a raster stack of abundance maps.
  # If the ID of the Replicate ("rep") is not provided, it will return the mean abundance over all replicates.
  
  if (!is.null(rep)){
    pop_wide <- reshape(subset(pop_df,Rep==rep)[,c('Year','x','y','NInd')], timevar='Year', v.names=c('NInd'), idvar=c('x','y'), direction='wide')
    r_years <- rasterFromXYZ(pop_wide)
    
    if (!is.null(mask)){
      r_years <- extend(r_years, mask)
      values(r_years)[is.na(values(r_years))] <- 0
      r_years <- mask(r_years, mask)
    }
    
  } else {
    pop_wide <- lapply(unique(pop_df$Year),FUN=function(year){reshape(subset(pop_df,Year==year)[,c('Rep','x','y','NInd')], timevar='Rep', v.names=c('NInd'), idvar=c('x','y'), direction='wide')})
    r_years <- stack(sapply(pop_wide, FUN=function(i){mean(extend(rasterFromXYZ(i),ext))}))
    names(r_years) <- paste0('mean.NInd.',unique(pop_df$Year))
    
    if (!is.null(mask)){
      r_years <- extend(r_years, mask)
      values(r_years)[is.na(values(r_years))] <- 0
      r_years <- mask(r_years, mask)
    }
  }
  return(r_years)
}


par(mfrow=c(1,1))

ext <-  c(min(pop_df$x)-100000,max(pop_df$x)+100000,min(pop_df$y)-100000,max(pop_df$y)+100000)

spain.boundry<- raster::shapefile("extras/Spain_mainland_map.shp")

# get the abundance maps with rep no-0
r_years_rep0 <- stack_pop(pop_df, ext, rep=2)
names(r_years_rep0)

result <- r_years_rep0[['NInd.10']]
crs(result) <- "EPSG:2062" 
#writeRaster(result,"Output_Maps/Spain_map_mink3cl_rep0_2020_R_1,2", overwrite=TRUE)

plot(spain.boundry)
plot(result,col=rev(brewer.pal(n = 5, name = "Reds")), axes=F, add =T )
plot(rasterToPolygons(SpDist, dissolve=F), add=T)



# get the abundance maps with mean of all repetition
r_years <- stack_pop(pop_df, ext)
names(r_years)
result <- r_years[['mean.NInd.20']]
crs(result) <- "EPSG:2062" 
#writeRaster(result,"Output_Maps/Spain_map_mink3cl_2020_R_1,2", overwrite=TRUE)

# plot the final mean abundance map and save it.
plot(spain.boundry)
plot(result,col=rev(brewer.pal(n = 5, name = "Reds")), axes=F, add =T )
plot(rasterToPolygons(SpDist, dissolve=F), add=T)




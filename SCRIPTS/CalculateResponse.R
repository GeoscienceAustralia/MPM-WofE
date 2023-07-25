#WofE - Calculate Weights Tool
#@author Arianne Ford, Geoscience Australia
#Based on ArcSDM Weights of Evidence toolbox

#get required packages for processing
require(tcltk)
require(rgdal)
require(terra)
require(vita)
require(caret)
require(sp)
require(dplyr)

#set working directory
workingdir <- tcltk::tk_choose.dir(caption="Select working directory for analysis")
setwd(workingdir)

#import shapefile
shapeFileFilter <- matrix(c("Shape Files", "*.shp"), ncol=2, byrow=T)
depositsPath <- tcltk::tk_choose.files(caption="Select mineral deposit training data", filters=shapeFileFilter, multi=F)

if (length(depositsPath) == 0 || !file.exists(depositsPath))
  stop("Could not locate mineral deposits shapefile...")

#read shapefile into spatVect
deposits <- terra::vect(depositsPath)

#get number of input maps for model and missing data value from user
askforInput <- function(inputMapNo){
  if (inputMapNo == "getMapNo"){
    getParam <- readline(prompt="How many input maps does your model have: ")
  }
  else if (inputMapNo == "missData"){
    getParam <- readline(prompt="Enter missing data value: ")
  }
  return(getParam)
}

#get number of input maps for model and missing data value from user
mapNo <- askforInput("getMapNo")
missData <- askforInput("missData")
mapNo <- as.numeric(mapNo)
missData <- as.numeric(missData)

for (mapNum in 1:mapNo){
  #iteratively import rasters and corresponding weights tables to make mineral potential map
  #import raster
  rasterFilter <- matrix(c("TIFF Image", "*.tif"), ncol=2, byrow=T)
  rasterPath <- tcltk::tk_choose.files(caption="Select integer raster for input to mineral potential map", filters=rasterFilter, multi=F)
  
  if (length(rasterPath) == 0 || !file.exists(rasterPath))
    stop("Could not locate raster file...")
 
  #read the raster into spatRast
  rastStack <- terra::rast(rasterPath)

  #make frequency table for raster for join to weights table
  freqTable <- freq(rastStack)
  
  #import corresponding weights table from calculate weights tool
  csvFilter <- matrix(c("CSV Weights Tables", "*.csv"), ncol=2, byrow=T)
  csvPath <- tcltk::tk_choose.files(caption="Select corresponding weights table to use for mineral potential map", filters=csvFilter, multi=F)
  
  #read csv into dataframe
  getCSV <- read.csv(csvPath)
  
  print("Processing...")
  
  #join csv to raster frequency table and read into matrix for analysis
  csvRast <- merge(freqTable, getCSV, by.x = "value", by.y = "CLASS", all.x = TRUE)
  data.matrix(csvRast)
  
  #get relevant values for calculations
  s = (sum(csvRast[,7]))
  ds = nrow(deposits)
  prior_prob = ds/s
  prior_logit = log(prior_prob/(1-prior_prob))
  wstd_const = 1.0/ds
  
  #calculate weighted raster for input to posterior logit calculation
  print(paste("Processing...Making weight raster for input map "), mapNum)
  for (rowNum in 1:nrow(freqTable)){
    if(rowNum==1){
      calcW <- terra::ifel(rastStack==csvRast[rowNum,1], (csvRast[rowNum,17]), (csvRast[rowNum+1,17]))
    }
    else{
      calcW <- terra::ifel(rastStack==csvRast[rowNum,1], (csvRast[rowNum,17]), calcW)
    }
  }

  #update posterior logit raster
  print(paste("Processing...Updating posterior logit with map "), mapNum)
  if(mapNum == 1){
    postlogit <- calcW + prior_logit
  }
  else{
    postlogit <- postlogit + calcW
  }

  #calculate standard deviation raster for input to weight standard deviation calculation
  print(paste("Processing...Updating STD raster for input map "), mapNum)
  for (rowNum in 1:nrow(freqTable)){
    if(rowNum==1){
      calcWSTD <- terra::ifel(rastStack==csvRast[rowNum,1], csvRast[rowNum,18], csvRast[rowNum+1,18])
    }
    else{
      calcWSTD <- terra::ifel(rastStack==csvRast[rowNum,1], csvRast[rowNum,18], calcWSTD)
    }
  }
  
  #update weight standard deviation
  if(mapNum == 1){
    sumWSTD <- (wstd_const) + (calcWSTD^2)
  }
  else{
    sumWSTD <- sumWSTD + (calcWSTD^2)
  }
}

#calculate posterior probability, standard deviation of posterior probaility, and posterior probability confidence rasters
postodds <- exp(postlogit)
print("Making posterior probability raster")
postprob <- postodds/(1+postodds)

print("Making standard deviation raster")
pprb_std <- sqrt((postprob^2) * sumWSTD)

print("Making confidence raster")
pprb_conf <- postprob/pprb_std

#save output rasters 
rasterFilter <- paste('{"TIFF Image" {"*.tif"}} ', sep=" ")
PostProb_out <- tcltk::tkgetSaveFile(title="Save posterior probability raster as",
                                       filetypes=rasterFilter)
if (nchar(PostProb_out <- as.character(PostProb_out)) == 0)
  stop("No raster file selected for saving posterior probability raster.")

writeRaster(postprob, PostProb_out, overwrite = TRUE)

PostProbSTD_out <- tcltk::tkgetSaveFile(title="Save standard deviation raster as",
                                     filetypes=rasterFilter)
if (nchar(PostProbSTD_out <- as.character(PostProbSTD_out)) == 0)
  stop("No raster file selected for saving posterior probability STD raster.")

writeRaster(pprb_std, PostProbSTD_out, overwrite = TRUE)

PostProbConf_out <- tcltk::tkgetSaveFile(title="Save confidence raster as",
                                     filetypes=rasterFilter)
if (nchar(PostProbConf_out <- as.character(PostProbConf_out)) == 0)
  stop("No raster file selected for saving posterior probability confidence raster.")

writeRaster(pprb_conf, PostProbConf_out, overwrite = TRUE)
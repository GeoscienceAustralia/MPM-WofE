#WofE - Calculate Weights Tool
#@author Arianne Ford, Geoscience Australia
#Based on ArcSDM Weights of Evidence toolbox (https://github.com/gtkfi/ArcSDM)

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

#import raster
rasterFilter <- matrix(c("TIFF Image", "*.tif"), ncol=2, byrow=T)
rasterPath <- tcltk::tk_choose.files(caption="Select integer raster to test", filters=rasterFilter, multi=F)
base_name <- basename(rasterPath)
raster_filename <- tools::file_path_sans_ext(base_name)

if (length(rasterPath) == 0 || !file.exists(rasterPath))
  stop("Could not locate raster file...")

# Read the image into spatRast
unclassRast <- terra::rast(rasterPath)
NAflag(unclassRast) <- -1

if(is.int(unclassRast) == FALSE)
  stop("Not an integer raster...")

if(is.lonlat(unclassRast) == TRUE)
  stop("Please use a projected coordinate system in meters!")

#import shapefile
shapeFileFilter <- matrix(c("Shape Files", "*.shp"), ncol=2, byrow=T)
depositsPath <- tcltk::tk_choose.files(caption="Select mineral deposit training data", filters=shapeFileFilter, multi=F)

if (length(depositsPath) == 0 || !file.exists(depositsPath))
  stop("Could not locate mineral deposits shapefile...")

#read shapefile into spatVect
deposits <- terra::vect(depositsPath)

# if(terra::same.crs(deposits, unclassRast) == FALSE){
#   stop("Raster and deposits must be in the same coordinate system")
# }

#request input from user
askforInput <- function(inputType){
  #is calculation to be ascending, descending, or categorical?
  if(inputType == "calcType"){
    getParam <- readline(prompt="Enter calculation type ((a)scending, (d)escending, (c)categorical: ")
  }
  #what is the confidence level required for studC to be valid?
  else if (inputType == "confCont"){
    getParam <- readline(prompt="Enter confidence required for studentized contrast: ")
  }
  #what is the unit area?
  else if(inputType == "unitArea"){
    getParam <- readline(prompt="Enter unit area (km2): ")
  }
  #what is the missing data value? 
  #This is not the same as NA values outside the study area mask
  else if (inputType == "missData"){
    getParam <- readline(prompt="Enter missing data value: ")
  }
  return(getParam)
}

calcType <- askforInput("calcType")
confCont <- askforInput("confCont")
unitArea <- askforInput("unitArea")
missData <- askforInput("missData")

calcType <- as.character(calcType)
confCont <- as.numeric(confCont)
unitArea <- as.numeric(unitArea)
missData <- as.numeric(missData)

# Create frequency table of cell counts
freqTable <- freq(unclassRast)

# Use sample points to extract raster signatures
depSig <- terra::extract(unclassRast, deposits, fun=NULL, method="simple", cells=FALSE, xy=FALSE)
colnames(depSig)[2] <- "value"
valCount <- merge(freqTable, depSig, by.x = "value", by.y="value")
by_val <- valCount %>% group_by(value)
valCount2 <- summarise(by_val, n = n())

#generate updated frequency tables for calculations
freqTable2 <- merge(freqTable, valCount2, by.x = "value", by.y="value", all.x=TRUE)
freqTable2[is.na(freqTable2)] <- 0
freqTableD <- freqTable2[order(freqTable2[,1], decreasing=TRUE),]

# get size of cell in km for each raster class
class_km2 = (unclassRast@ptr[["res"]])/(1000)

# Make weights table
makeWts <- matrix(nrow=nrow(freqTable), ncol=15)
colnames(makeWts) = c("CLASS", "CLASS_AREA", "AREA_SQ_KM", "AREA_UNITS", "NO_POINTS", "WPLUS", "S_WPLUS", "WMINUS", "S_WMINUS", "CONTRAST", "S_CONTRAST", "STUDC", "GEN_CLASS", "WEIGHT", "W_STD")

#calculations if calculation type is "descending" from user
if (calcType == "d" | calcType == "D"){
  for (row in 1:nrow(makeWts)){
    makeWts[row, 1] = freqTableD[row, 1]
    #set map values and make area calculations
    if ((row == 1) & (makeWts[row, 1] != missData)){
      makeWts[row, 2] = freqTableD[row, 3] * class_km2[1] * class_km2[1]
      makeWts[row, 3] = freqTableD[row, 3] * class_km2[1] * class_km2[1]
      makeWts[row, 4] = makeWts[row, 3] * unitArea
      makeWts[row, 5] = freqTableD[row, 4]
    }
    if ((row > 1) & (makeWts[row,1] != missData)){
      makeWts[row, 2] = (freqTableD[row, 3] * class_km2[1] * class_km2[1])
      makeWts[row, 3] = (freqTableD[row, 3] * class_km2[1] * class_km2[1]) + makeWts[(row-1),3]
      makeWts[row, 4] = (makeWts[row, 3]) * unitArea
      makeWts[row, 5] = freqTableD[row, 4] + makeWts[(row-1), 5]
    }
    #set all values to 0 if map value is equal to missing data value from user
    if (makeWts[row, 1] == missData){
      makeWts[row, 2] = 0.0
      makeWts[row, 3] = 0.0
      makeWts[row, 4] = 0.0
      makeWts[row, 5] = 0.0
      makeWts[row, 6] = 0.0
      makeWts[row, 7] = 0.0
      makeWts[row, 8] = 0.0
      makeWts[row, 9] = 0.0
      makeWts[row, 10] = 0.0
      makeWts[row, 11] = 0.0
      makeWts[row, 12] = 0.0
      makeWts[row, 13] = missData
      makeWts[row, 14] = 0.0
      makeWts[row, 15] = 0.0
    }
  }
  
  #calculate statistics for weights table
  for (row in 1:nrow(makeWts)){
    s = max(makeWts[, 4])
    b = makeWts[row, 4]
    ds = nrow(deposits)
    db = makeWts[row, 5]
    
    #calculate prior probability
    priorprob = ds/s
    
    #error checking and patches for calculations to avoid divide by 0 errors
    if(db>ds){
      print("Input error: more than one training point per unit cell in study area!")
      return ()
    }
    if(db==0){ #no deposit accumulation
      makeWts[row, 6] = 0.0
      makeWts[row, 7] = 0.0
      makeWts[row, 8] = 0.0
      makeWts[row, 9] = 0.0
      makeWts[row, 10] = 0.0
      makeWts[row, 11] = 0.0
      makeWts[row, 12] = 0.0
    }
    if(db==ds){ #maximum accumulation
      db = db-0.99
    }
    if((s-b) <= (ds-db)){ #fix for db = maxTP
      b = s+db-ds-0.99
    } 
    if((b-db) <= 0.0){ #fix if area < unit size
      b = db+1
    }
    
    #calculate statistics for cumulative map value
    if (makeWts[row, 1] != missData & db!=0){
      #calculate W+
      pbd = db/ds
      pbdb = (b-db) / (s-ds)
      ls = pbd/pbdb
      wp = log(ls)
      makeWts[row, 6] = wp
      
      #calculate vp and sp
      vp = (1.0 / db) + (1.0 / (b-db))
      sp = sqrt(vp)
      makeWts[row, 7] = sp
      
      #calculate W-
      pbbd = (ds-db) / ds
      pbbdb = (s-b-ds+db) / (s-ds)
      ln1 = pbbd / pbbdb
      wm = log(ln1)
      makeWts[row, 8] = wm
      
      #calculate vm and sm
      vm = (1.0 / (ds-db)) + (1.0 / (s-b-ds+db))
      sm = sqrt(vm)
      makeWts[row, 9] = sm
      
      #calculate contrast
      c = wp - wm
      makeWts[row, 10] = c
      
      #calculate contrast stddev and studC
      sc = sqrt(vp+vm)
      makeWts[row, 11] = sc
      makeWts[row, 12] = c/sc
    }
    
    #reset if missing data value 
    if (makeWts[row, 1] == missData){
      makeWts[row, 6] = 0.0
      makeWts[row, 7] = 0.0
      makeWts[row, 8] = 0.0
      makeWts[row, 9] = 0.0
      makeWts[row, 10] = 0.0
      makeWts[row, 11] = 0.0
      makeWts[row, 12] = 0.0
    }
  }
  
  for (row in 1:nrow(makeWts)){  
    #determine genclass and set generalised weights for gen class
    maxConf = which.max(makeWts[, 12])
    maxCont = which.max(makeWts[, 10])
    if (makeWts[maxConf,12] < confCont){
      makeWts[, 13] = 999
      makeWts[, 14] = 0.0
      makeWts[, 15] = 0.0
    }
    else if ((makeWts[row, 1]<=max(makeWts[,1])) & (makeWts[row, 10]<=max(makeWts[,10])) & (makeWts[maxConf,12] >= confCont) & (row<=maxCont)){
      makeWts[row, 13] = 2
      makeWts[row, 14] = (makeWts[maxCont, 6])
      makeWts[row, 15] = (makeWts[maxCont, 7])
    }
    else if ((makeWts[row, 1]<=max(makeWts[,1])) & (makeWts[row, 10]<=max(makeWts[,10])) & (makeWts[maxConf,12] >= confCont) & (row> maxCont)){
      makeWts[row, 13] = 1
      makeWts[row, 14] = (makeWts[maxCont, 8])
      makeWts[row, 15] = (makeWts[maxCont, 9])
    }
    if (makeWts[row, 1] == missData){
      makeWts[row, 13] = missData
      makeWts[row, 14] = 0.0
      makeWts[row, 15] = 0.0
    }
  }
}

#calculations if calculation type is "ascending" from user
if (calcType == "a" | calcType == "A" ){
  for (row in 1:nrow(makeWts)){
    makeWts[row, 1] = freqTable[row, 2]
    #set map values and make area calculations
    if ((row == 1) & (makeWts[row,1]!= missData)){
      makeWts[row, 2] = freqTable[row, 3] * class_km2[1] * class_km2[1]
      makeWts[row, 3] = freqTable[row, 3] * class_km2[1] * class_km2[1]
      makeWts[row, 4] = makeWts[row, 3] * unitArea
      makeWts[row, 5] = freqTable2[row, 4]
    }
    if ((row > 1) & (makeWts[row,1] != missData)){
      makeWts[row, 2] = (freqTable[row, 3] * class_km2[1] * class_km2[1])
      makeWts[row, 3] = (freqTable[row, 3] * class_km2[1] * class_km2[1]) + makeWts[(row-1),3]
      makeWts[row, 4] = (makeWts[row, 3] * unitArea) 
      makeWts[row, 5] = freqTable2[row, 4] + makeWts[(row-1), 5]
    }
    #set all values to 0 if map value is equal to missing data value from user
    if (makeWts[row, 1] == missData){
      makeWts[row, 2] = 0.0
      makeWts[row, 3] = 0.0
      makeWts[row, 4] = 0.0
      makeWts[row, 5] = 0.0
      makeWts[row, 6] = 0.0
      makeWts[row, 7] = 0.0
      makeWts[row, 8] = 0.0
      makeWts[row, 9] = 0.0
      makeWts[row, 10] = 0.0
      makeWts[row, 11] = 0.0
      makeWts[row, 12] = 0.0
      makeWts[row, 13] = missData
      makeWts[row, 14] = 0.0
      makeWts[row, 15] = 0.0
    }
  }
  
  #calculate statistics for weights table
  for (row in 1:nrow(makeWts)){
    s = max(makeWts[, 4])
    b = makeWts[row, 4]
    ds = nrow(deposits)
    db = makeWts[row, 5]
    
    #calculate prior probability
    priorprob = ds/s

    #error checking and patches for calculations to avoid divide by 0 errors
    if(db>ds){
      print("Input error: more than one training point per unit cell in study area!")
      return ()
    }
    if(db==0){ #no accumulation
      makeWts[row, 6] = 0.0
      makeWts[row, 7] = 0.0
      makeWts[row, 8] = 0.0
      makeWts[row, 9] = 0.0
      makeWts[row, 10] = 0.0
      makeWts[row, 11] = 0.0
      makeWts[row, 12] = 0.0
    }
    if(db==ds){ #maximum accumulation
      db = db-0.99
    }
    if((s-b) <= (ds-db)){ #fix for db = maxTP
      b = s+db-ds-0.99
      } 
    if((b-db) <= 0.0){ #fix if area < unit size
      b = db+1
    }
    
    #calculate statistics for cumulative map value
    if (makeWts[row, 1] != missData & db!=0){
      #calculate W+
      pbd = db/ds
      pbdb = (b-db) / (s-ds)
      ls = pbd/pbdb
      wp = log(ls)
      makeWts[row, 6] = wp
      
      #calculate vp and sp
      vp = (1.0 / db) + (1.0 / (b-db))
      sp = sqrt(vp)
      makeWts[row, 7] = sp
      
      #calculate W-
      pbbd = (ds-db) / ds
      pbbdb = (s-b-ds+db) / (s-ds)
      ln1 = pbbd / pbbdb
      wm = log(ln1)
      makeWts[row, 8] = wm
      
      #calculate vm and sm
      vm = (1.0 / (ds-db)) + (1.0 / (s-b-ds+db))
      sm = sqrt(vm)
      makeWts[row, 9] = sm
      
      #calculate contrast
      c = wp - wm
      makeWts[row, 10] = c
      
      #calculate contrast stddev and studC
      sc = sqrt(vp+vm)
      makeWts[row, 11] = sc
      makeWts[row, 12] = c/sc
    }
    
    #reset if missing data value 
    if (makeWts[row, 1] == missData){
      makeWts[row, 6] = 0.0
      makeWts[row, 7] = 0.0
      makeWts[row, 8] = 0.0
      makeWts[row, 9] = 0.0
      makeWts[row, 10] = 0.0
      makeWts[row, 11] = 0.0
      makeWts[row, 12] = 0.0
    }
  }
  
  for (row in 1:nrow(makeWts)){ 
    #determine genclass and set generalised weights for gen class
    maxConf = which.max(makeWts[, 12])
    maxCont = which.max(makeWts[, 10])
    if (makeWts[maxConf,12] < confCont){
      makeWts[, 13] = 999
      makeWts[, 14] = 0.0
      makeWts[, 15] = 0.0
    }
    else if ((makeWts[row, 1]<=max(makeWts[,1])) & (makeWts[row, 10]<=max(makeWts[,10])) & (makeWts[maxConf, 12] >= confCont) & (row<=maxCont)){
      makeWts[row, 13] = 2
      makeWts[row, 14] = (makeWts[maxCont, 6])
      makeWts[row, 15] = (makeWts[maxCont, 7])
    }
    else if ((makeWts[row, 1]<=max(makeWts[,1])) & (makeWts[row, 10]<=max(makeWts[,10])) & (makeWts[maxConf,12] >= confCont) & (row> maxCont)){
      makeWts[row, 13] = 1
      makeWts[row, 14] = (makeWts[maxCont, 8])
      makeWts[row, 15] = (makeWts[maxCont, 9])
    }
    if (makeWts[row, 1] == missData){
      makeWts[row, 13] = missData
      makeWts[row, 14] = 0.0
      makeWts[row, 15] = 0.0
    }
  }
}

#calculations if calculation type is "categorical" from user
if (calcType == "c" | calcType == "C"){
  for (row in 1:nrow(makeWts)){
    makeWts[row, 1] = freqTable2[row, 1]
    #set map values and make area calculations
    if ((makeWts[row,1]!= missData)){
      makeWts[row, 2] = freqTable2[row, 3] * class_km2[1] * class_km2[1]
      makeWts[row, 3] = freqTable2[row, 3] * class_km2[1] * class_km2[1]
      makeWts[row, 4] = freqTable2[row, 3] * class_km2[1] * class_km2[1] * unitArea
      makeWts[row, 5] = freqTable2[row, 4]
    }
    
    #set all values to 0 if map value is equal to missing data value from user
    if (makeWts[row, 1] == missData){
      makeWts[row, 2] = 0.0
      makeWts[row, 3] = 0.0
      makeWts[row, 4] = 0.0
      makeWts[row, 5] = 0.0
      makeWts[row, 6] = 0.0
      makeWts[row, 7] = 0.0
      makeWts[row, 8] = 0.0
      makeWts[row, 9] = 0.0
      makeWts[row, 10] = 0.0
      makeWts[row, 11] = 0.0
      makeWts[row, 12] = 0.0
      makeWts[row, 13] = missData
      makeWts[row, 14] = 0.0
      makeWts[row, 15] = 0.0
    }
  }
  
  #calculate statistics for weights table
  for (row in 1:nrow(makeWts)){
    s = as.numeric(sum(makeWts[, 4]))
    b = as.numeric(makeWts[row, 4])
    ds = nrow(deposits)
    db = as.numeric(makeWts[row, 5])
    
    #calculate prior probability
    priorprob = ds/s

    #error checking and patches for calculations to avoid divide by 0 errors
    if(db>ds){
      print("Input error: more than one training point per unit cell in study area!")
      return ()
    }
    if(db==0){ #no accumulation
      makeWts[row, 6] = 0.0
      makeWts[row, 7] = 0.0
      makeWts[row, 8] = 0.0
      makeWts[row, 9] = 0.0
      makeWts[row, 10] = 0.0
      makeWts[row, 11] = 0.0
      makeWts[row, 12] = 0.0
    }
    if(db==ds){ #maximum accumulation
      db = db-0.99
    }
    if(db==0.001){ #fix for old ArcSDM problem, can possibly now be removed?
      db = ds
      db = db-0.99
    }
    if((s-b) <= (ds-db)){ #fix for db = maxTP
      b = s+db-ds-0.99
    } 
    if((b-db) <= 0.0){ #fix if area < unit size
      b = db+1
    }
    
    #calculate statistics for individual (non-cumulative) categorical map values
    if (makeWts[row, 1] != missData & db!=0){
      #calculate W+
      pbd = db/ds
      pbdb = (b-db) / (s-ds)
      ls = pbd/pbdb
      wp = log(ls)
      makeWts[row, 6] = wp
      
      #calculate vp and sp
      vp = (1.0 / db) + (1.0 / (b-db))
      sp = sqrt(vp)
      makeWts[row, 7] = sp
      
      #calculate W-
      pbbd = (ds-db) / ds
      pbbdb = (s-b-ds+db) / (s-ds)
      ln1 = pbbd / pbbdb
      wm = log(ln1)
      makeWts[row, 8] = wm
      
      #calculate vm and sm
      vm = (1.0 / (ds-db)) + (1.0 / (s-b-ds+db))
      sm = sqrt(vm)
      makeWts[row, 9] = sm
      
      #calculate contrast
      c = wp - wm
      makeWts[row, 10] = c
      
      #calculate contrast stddev and studC
      sc = sqrt(vp+vm)
      makeWts[row, 11] = sc
      makeWts[row, 12] = c/sc
    }
    
    #reset if missing data value 
    if (makeWts[row, 1] == missData){
      makeWts[row, 6] = 0.0
      makeWts[row, 7] = 0.0
      makeWts[row, 8] = 0.0
      makeWts[row, 9] = 0.0
      makeWts[row, 10] = 0.0
      makeWts[row, 11] = 0.0
      makeWts[row, 12] = 0.0
    }
  }
  
  for (row in 1:nrow(makeWts)){
    #determine genclass and set generalised weights for gen class
    maxConf = which.max(makeWts[, 12])
    if (makeWts[maxConf,12] < confCont){
      makeWts[, 13] = 999
      makeWts[, 14] = 0.0
      makeWts[, 15] = 0.0
    }
    else if ((makeWts[row, 12] >= confCont)){
      makeWts[row, 13] = makeWts[row, 1]
      makeWts[row, 14] = makeWts[row, 6]
      makeWts[row, 15] = makeWts[row, 7]
    }
    else if ((makeWts[row, 12] < confCont) & makeWts[maxConf,12] >= confCont){
      makeWts[row, 13] = makeWts[row, 1]
      makeWts[row, 14] = makeWts[row, 6]
      makeWts[row, 15] = makeWts[row, 7]
    }
    if (makeWts[row, 1] == missData){
      makeWts[row, 13] = missData
      makeWts[row, 14] = 0.0
      makeWts[row, 15] = 0.0
    }
  }
}

print(paste("Prior probability = ", priorprob))

#save weights tables as csv with appropriate suffix depending on calculation type
if (calcType == "d" | calcType == "D"){
  CD_Type <- paste(raster_filename, "CD", sep="_")
  writeWts <- paste(CD_Type, "csv", sep=".")
  write.csv(makeWts, file = writeWts)
}

if (calcType == "a" | calcType == "A"){
  CA_Type <- paste(raster_filename, "CA", sep="_")
  writeWts <- paste(CA_Type, "csv", sep=".")
  write.csv(makeWts, file = writeWts)
}

if (calcType == "c" | calcType == "C"){
  CT_Type <- paste(raster_filename, "CT", sep="_")
  writeWts <- paste(CT_Type, "csv", sep=".")
  write.csv(makeWts, file = writeWts)
}

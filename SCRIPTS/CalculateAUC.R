#WofE - Area Under the Curve Tool
#@author Arianne Ford, Geoscience Australia
#Based on ArcSDM Weights of Evidence (Area-Frequency) tool (https://github.com/gtkfi/ArcSDM)

#get required packages for processing
library(modEvA)
library(terra)
library(rgdal)
library(sp)
library(pROC)

#set working directory
workingdir <- tcltk::tk_choose.dir(caption="Select working directory for analysis")
setwd(workingdir)

#import raster
rasterFilter <- matrix(c("TIFF Image", "*.tif"), ncol=2, byrow=T)
rasterPath <- tcltk::tk_choose.files(caption="Select raster to test", filters=rasterFilter, multi=F)

if (length(rasterPath) == 0 || !file.exists(rasterPath))
  stop("Could not locate raster file...")

# Read the image into dataframe for analysis
minpot <- rast(rasterPath)

#import shapefile
shapeFileFilter <- matrix(c("Shape Files", "*.shp"), ncol=2, byrow=T)
depositsPath <- tcltk::tk_choose.files(caption="Select mineral deposit training data", filters=shapeFileFilter, multi=F)

if (length(depositsPath) == 0 || !file.exists(depositsPath))
  stop("Could not locate mineral deposits shapefile...")

#read shapefile into dataframe for analysis
deposits <- readOGR(dirname(depositsPath), tools::file_path_sans_ext(basename(depositsPath)), integer64="no.loss", GDAL1_integer64_policy=FALSE)

#get coordinates of deposits for AUC calculations
depcoords <- coordinates(deposits)

#calculate AUC and plot graph
minpotAUC <- AUC(obs = depcoords, pred = minpot, main = "ROC curve")

#export AUC statistics to csv
csvFilter <- paste('{"CSV file" {"*.csv"}} ', sep=" ")
AUC_out <- tcltk::tkgetSaveFile(title="Save AUC statistics as",
                                     filetypes=csvFilter)
if (nchar(AUC_out <- as.character(AUC_out)) == 0)
  stop("No raster file selected for saving AUC statistics.")

write.csv(minpotAUC, AUC_out)

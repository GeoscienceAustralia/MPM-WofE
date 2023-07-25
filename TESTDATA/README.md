# Test Data for WofE R scripts

Please store files in a folder path that does not contain any spaces in the pathname!

## Test Maps
### Calculate Weights
Four maps have been derived over the Cloncurry IOCG district based on data in: 
* Stewart, A.J., Liu, S.F., Highet, L.M., Woods, M., Czarnota, K., Bonnardot, M., Brown, C., Clark, A., Connors, K. 2020. Solid Geology of the North Australian Craton, 1:1 000 000 scale, 1st edition (2020). Geoscience Australia, Canberra. https://dx.doi.org/10.26186/135277

The maps are Geotiff files with an integer data type (required):
* Dist2Breccia.tif - Distance to brecciated units queried from the solid geology. Should be tested using ascending calculation type because low values are favourable.
* Dist2Fault.tif - Distance to faults (no subset applied). Should be tested using ascending calculation type because low values are favourable.
* FaultDensity.tif - Density of faults (no subset applied) using a 25km search radius and converted to integer grid using 10 natural breaks. Should be tested using descending calculation type because high values are favourable.
* StratNo.tif - Solid geology reclassified by the StratNo (stratigraphic unit number) field. Should be tested using categorical calculation type because the StratNo values are unordered in terms of favourability.

### Calculate Response
The binary maps in the Calculate Response folder use the following thresholds:
* Dist2Breccia.tif - 5km distance buffer
* Dist2Fault.tif - 2km distance buffer
* Fault Density.tif - Class 10-6 (0.371 - 0.186 faults/km2)
* StratNo.tif - StratNo = 76890 or 42515 or 31441 or 26091 or 12754 or 11361 or 7245 ()

The weights tables have been regenerated for the binary maps using the categorical calculation type, because the binary values are now unordered.

### Calculate AUC
The map in the Calculate AUC folder is the posterior probability (mineral potential) map generated from the Calculate Response script using the 4 binary input maps described above. For the purposes of the test data, the mineral deposits and occurrences have not been subset into training and validation sets.

## Training Data
The IOCG mineral deposits and occurrences used for testing and validating the maps have been derived from the following datasets:
* Huston, D., Doublier, M., Downes, P.M. 2021. Geological setting, age and endowment of major Australian mineral deposits - a compilation. Geoscience Australia, Canberra. http://dx.doi.org/10.11636/Record.2021.020
* Kucka, C., Senior, A., Britt, A. 2022. Mineral Occurrences: Forgotten discoveries providing new leads for mineral supply. Geoscience Australia, Canberra. https://dx.doi.org/10.26186/146983

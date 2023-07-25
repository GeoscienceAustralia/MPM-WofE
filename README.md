# Mineral Potential Mapping - Weights of Evidence

Weights of evidence (WofE) is a data-driven method for mineral potential mapping that is based on Bayesian statistics. It evaluates the relationship between a set of known mineral deposits/occurrences and an individual map and compares the observed correlation to the expected correlation. 

Three R scripts have been developed as part of Geoscience Australia's Exploring for the Future program. Geoscience Australia’s Exploring for the Future program provides precompetitive information to inform decision-making by government, community, and industry on the sustainable development of Australia's mineral, energy, and groundwater resources. By gathering, analysing, and interpreting new and existing precompetitive geoscience data and knowledge, we are building a national picture of Australia’s geology and resource potential. This leads to a strong economy, resilient society and sustainable environment for the benefit of all Australians. This includes supporting Australia’s transition to net zero emissions, strong, sustainable resources and agriculture sectors, and economic opportunities and social benefits for Australia’s regional and remote communities. The Exploring for the Future program, which commenced in 2016, is an eight year, $225m investment by the Australian Government.

Based on the ArcSDM toolbox for ArcGIS Pro (https://github.com/gtkfi/ArcSDM), three R scripts have been developed that replicate the Calculate Weights, Calculate Response, and Area-Frequency tools within the Weights of Evidence toolbox. The R scripts have been developed in order to remove the dependency on ArcGIS and the Spatial Analyst toolbox.

Calculate Weights evaluates the Contrast and Studentized contrast values for each map value and indicates where the optimal threshold occurs using the GEN_CLASS field. It is up to the user to determine if that threshold is geologically meaningful, and manually change it if required.

Calculate Response generates the mineral potential map for a user specified number of input maps and their corresponding weights tables from the Calculate Weights tool. The script also produces a corresponding standard deviation and confidence map.

Calculate AUC evaluates the area under the curve for an individual input map or output mineral potential map. It outputs an Receiver-Operating Characteristic curve and a table containing relevant statistics for additional analysis.

# Dependencies

The scripts have been developed in R Studio 2022.02.2 (Build 485) using the following packages and versions:
* tcltk (1.2-11)
* rgdal (1.6-2)
* terra (1.5-34)
* vita (1.0.0)
* caret (6.0-93)
* sp (1.5-0)
* dplyr (1.0.9)
* modEvA (3.82)
* pROC (1.18.0)

# Running

## Calculate Weights
* Calculate Weights requires a Geotiff file (must be Integer type) that represents the map to be tested and a shapefile containing the known mineral deposits/occurrrences. Both the Geotiff and shapefile must be in the same coordinate system. The script assumes a projected coordinate system in meters.
* The script requests information from the user including:
	* Calculation type - use ascending where low map values are favourable, use descending when high map values are favourable, use categorical for unordered map data such as lithology.
	* Confidence threshold for Studentized Contrast - minimum confidence required for results to be considered statistically valid (values between 0.5 and 2 are recommended). 
	* Unit Area in km2 - what is the estimated footprint size of the deposits (1km2 is recommended for regional-scale analysis)?
	* Missing Data Value - what is the missing data value (to account for missing data in the map)? This is not the same as NULL or N/A values that fall outside the study area mask.
* The Calculate Weights script outputs a weights table with the relevant statistics for user review.

## Calculate Response
* Requires a user specified number of Geotiff files (must be Integer type) that represent the maps to be integrated to produce the mineral potential map, along with each map's corresponding weights table from the Calculate Weights script, and a shapefile containing the known mineral deposits/occurrrences. Note that the Geotiff files and shapefile must be in the same coordinate system. The script assumes a projected coordinate system in meters.
* It is strongly recommended that users convert each input map to binary prior to input to the Calculate Response tool in order to speed up processing. The binary map threshold should be determined using the GEN_CLASS field in the Calculate Weights output table (or modified accordingly if the statistical threshold determined by the script does not make geological sense). The Calculate Weights tool should be re-run on the binary map prior to running Calculate Response.
* The script requests information from the user including:
	* Number of input maps - how many input maps are to be included in the mineral potential map?
	* Missing Data Value - what is the missing data value (to account for missing data in the map). This is not the same as NULL or N/A values that fall outside the study area mask.
* Calculate Response outputs 3 maps:
	* Posterior probability - the mineral potential map 
	* Standard deviation of posterior probability - standard deviation of the mineral potential map
	* Confidence in posterior probability - posterior probability/standard deviation of posterior probability

## Calculate AUC
* Calculate AUC requires an input Geotiff file (on a 0-1 probability scale) that represents the map to be tested and a shapefile containing the known mineral deposits/occurrences. Both the Geotiff and shapefile must be in the same coordinate system.
* The script outputs the Receiver-Operating Characteristic curve with the corresponding Area Under the Curve value. A csv file is also produced that contains statistics relevant to the AUC calculations.

# ELC10

**The methodology in this repository is currently under review at an academic journal. The link to the paper will be posted in this README when available. Please refer to the manuscript for details on methods**

The preprint article can be accessed here: https://arxiv.org/abs/2104.10922

This is workflow to classify land cover over Europe at 10 m resolution produced with Sentinel optical and radar satellite imagery. The classification model is trained on land cover reference data form the LUCAS (Land Use/Cover Area frame Survey) dataset. The map represents conditions in 2018.

The workflow spans two programming platforms, namely Google Earth Engine (GEE) JavaScript API and RStudio. The R code is merely to finalize and evaluate the classification model so that you know what you are getting when you make final inference in GEE. Careful attention needs to be paid to the movement of data between the two. The data to be uploaded to GEE are located int he ./DATA/For_GEE folder and the data to be used in the R script are found in ./DATA/From_GEE folder. 

Between each step in the workflow, please rename file paths in the JavaScript code which are currently hardcoded to a demo GEE Asset. You need to generate data within your own GEE Assets. Please also create Google Drive folders for the data to export into - these are described in-line in the code. The R code should not need to be edited.

The workflow is as follows:

- Upload the following files to your GEE Assets folder:
  - "./DATA/For_GEE/EEA_grid.shp" - this is a 100x100km grid used to break up processing into smaller chunks
  - "./DATA/For_GEE/ELC10_training_feats.csv" - this contains Lat and Lon coordinates for LUCAS reference points which have land cover labels. Please see manuscript for details on selection of these points.
- Run the following scripts in GEE:
  - script 'auxiliary_extarct.js' - this extracts climate, nighttime lights, and terrain variables
  - 'sentinel_2_extract.js' - this extracts Sentinel-2 spectral-temporal metrics
  - 'sentine_1_extract.js' - this extarcts Sentinel-1 spectral-temporal metrics
*wait for exports tasks to complete running*
- Download exported data from your Google Drive to the Rproject folder './DATA/From_GEE/' - pre-baked data area already in there. But replace if necessary.
- Run the R script 'RF_evaluation.R' - this trains and evaluates a Random Forest classification model and exports a final training dataset for upload to GEE.
- Upload the file "./DATA/For_GEE/LUCAS_combined_reference_final.csv" file to GEE Assets folder
- Run the 'RF_predict.js' script in GEE - this trains the GEE RF model and makes predictions over entire Europe.
*wait for RF prediction export tasks to complete running (3/4 days)*
- Run the 'Output_visualize.js' - this simply visualizes the classified land cover image.



*As it stands, this workflow is coded for efficiency because it is not performing speckle filtering on the Sentinel-1 data. For final implementation, and optimal accuracy, please perform speckle filtering*

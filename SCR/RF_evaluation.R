#### Setup environment ------------------------------------------------------------------------
# Import libraries (install if not already installed)
library(tidyverse)
library(randomForest)
library(caret)
library(sf)

# Set plotting theme for ggplot
theme_set(theme_bw()+ 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
            theme(strip.background =element_rect(fill="white")))

# Function to read in multiple files from directory
readMultiFiles <- function(directory){
  
  files <- list.files(directory, pattern='*.csv', full.names=TRUE)
  raw <- files %>% 
    map_df(~read_csv(.))
  return (raw)
  
}

#### Import data ------------------------------------------------------------------------------
# LUCAS polygon & point land cover reference points
ref_points <- read_csv('./DATA/For_GEE/LUCAS_combined_reference_final.csv') %>%
  dplyr::select(-votes)

# Sentinel 1 SAR data from GEE
s1 <- readMultiFiles('./DATA/From_GEE/Sentinel/S1_data/') %>%
  dplyr::select(-'system:index', -'.geo', -LC, -LC_num, -votes, -source, -X1)

# Sentinel 2 optical data from GEE
s2 <- readMultiFiles('./DATA/From_GEE/Sentinel/S2_data/') %>%
  dplyr::select(-'system:index', -'.geo', -LC, -LC_num, -votes, -source, -X1)

# Auxiliary data from GEE
auxiliary <- read_csv('./DATA/From_GEE/Auxiliary/auxiliary_data.csv') %>%
  dplyr::select(-'system:index', -'.geo', -LC, -LC_num, -votes, -source)

# Merge into one data frame
trainingFeats <- ref_points %>%
  left_join(s1, by="POINT_ID")%>%
  left_join(s2, by="POINT_ID")%>%
  left_join(auxiliary, by="POINT_ID") %>%
  mutate(LC = factor(LC, levels = c('Artificial land','Cropland', 'Woodland', 'Grassland',
                                    'Shrubland', 'Bare land', 'Water', 'Wetland'))) %>%
  drop_na()

#### Data checking and explore ------------------------------------------------------------
# See distribution of reference data
  #LUCAS_point are from the master LUCAS database
  #LUCAS_polygon are from the copernicus module
ref_points %>%
  ggplot(aes(x = LC, fill=source)) +
  geom_bar(stat='count') +
  coord_flip()

# Number of training features by data source
ref_points %>%
  group_by(source) %>%
  summarise(n=n())
nrow(ref_points)

# Check the names of all variables
names(trainingFeats)

#### Identify optimal number of predictors -------------------------------------------------
# The purpose of this is to identify a reduced set of predictors so that the 
# RF inference in GEE is efficient and can run faster

# List non-predictor variables
nonPredVars <- c('POINT_ID', 'LC', 'LC_num', 'source', 'lat', 'lon')

# List predictor variables
selectVars <- names(trainingFeats)[!names(trainingFeats) %in% nonPredVars]
selectVars
length(selectVars)

# Define Recursive Feature Selection (RFE) control parameters
control <- rfeControl(functions=rfFuncs, method="cv", number=5, verbose=TRUE)

# Define subsets for the function to iterate over
subsets <- c(2,5,10,15,20, 30, 40)

# run the RFE algorithm
results <- rfe(trainingFeats[,selectVars], trainingFeats$LC, sizes=subsets, rfeControl=control)

# summarize the results
print(results)

# Get the optimal numbe of variables
plot(results, type=c("g", "o")) # 15 variables is optimal


#### Identify top predictors to keep -------------------------------------------------------
# Now we rank the most important predictor variables and choose the top 15 for final model

# Create training dataset for RF model
training <-  trainingFeats %>%
  dplyr::select(all_of(c(selectVars, 'LC'))) %>%
  mutate(LC = factor(LC)) 

# Empty data frame for importance scores
import <- data.frame()
iteration <- seq(1,10, 1)

# Run for loop to generate importance scores
for (i in iteration){
  print(paste0('Iteration number: ', i))
  
  model.rf <-  randomForest(LC~.,
                            data=training, 
                            ntree=100,
                            importance=T, 
                            do.trace = 20)
  
  newImp <- as.data.frame(importance(model.rf))
  newImp$var <- rownames(newImp)
  newImp$iteration <- i
  import <- bind_rows(import, newImp)
}

import <- import %>% as_tibble()

# Calculate varibale importance as the average between normalized Gini and accuracy
topVars <- import %>%
  mutate(meanAcc = (BBmisc::normalize(MeanDecreaseAccuracy,  method = "range") + 
                      BBmisc::normalize(MeanDecreaseGini,  method = "range"))/2) %>%
  group_by(var) %>%
  summarise(MeanDecreaseGini = mean(MeanDecreaseGini),
            MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy),
            imp = mean(meanAcc), 
            sd = sd(meanAcc)/sqrt(n())) %>%
  arrange(desc(imp))

# Top 15 variables
topVars$var[1:15]

# PLot them
topVars %>%
  filter(var %in% topVars$var[1:15]) %>%
  ggplot(aes(x=reorder(var, imp), y=imp*100)) +
  geom_bar(stat='identity') + coord_flip() +
  ylab("Importance (%)") +
  xlab('')



#### Run final Random Forest -------------------------------------------------------------------
# Define final selected variables
#selectVarsFinal <- topVars$var[1:15]

# Hard coding in the top variables for illustrative purposes
  # These variables were derived from the final pre-processing workflow defined in the manuscript
selectVarsFinal <- c("ndvi_p25",
                     "green_median",
                     "temp" ,
                     "light",
                     "nbr_stdDev", 
                     "asc_vh_median" ,
                     "desc_vh_median",
                     "swir2_median" ,
                     "elevation" ,
                     "desc_dpol_median",
                     "R1_median",
                     "precip" ,
                     "asc_dpol_median",
                     "temp_stDev",
                     "asc_vv_stDev")

# Create final training dataset
trainingFinal <- trainingFeats %>%
  dplyr::select(all_of(c(selectVarsFinal, 'LC'))) %>%
  mutate(LC = factor(LC)) 

# Run RF model (ntree and default mtry have been tuned based on previous scripts)
model.rf <- randomForest(LC~.,
                         data=trainingFinal, 
                         ntree=100,
                         importance=T, 
                         do.trace = 20)
print(model.rf)
plot(model.rf)
varImpPlot(model.rf, sort=T)

# Confusion matrix
cm <- confusionMatrix(model.rf$predicted, trainingFinal$LC)
cm

# Export final training dataset for upload to GEE
trainingFeatsExport <-trainingFeats %>%
  dplyr::select(all_of(c( 'LC', 'LC_num', selectVarsFinal))) %>%
  mutate(LC = factor(LC))
trainingFeatsExport %>%
  write_csv('./DATA/For_GEE/ELC10_training_feats.csv')
 
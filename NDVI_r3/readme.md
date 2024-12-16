# NDVI data analysis using NNGP and other state-of-the-art methods

This folder provides scripts for this analysis. 
While we do not include full MCMC runs because they are too heavy, we put Gelman-Rubin-Brooks plots and trace plots in order to show that our chains have a sound behavior. 
They are located at */res*, under the *.pdf* format. 

## Run the model

If you want to replicate the NDVI experiment or you just want to harm your CPU for no particular reason, download the whole folder and follow those instructions. 
The data is in *data_cleaned_small_expanded.Rdata*. 
Then, open *run.R* and execute the script. 
The first part of the script will do a train-test partition. It will output the image *train_test_split.pdf* and several files corresponding to the train and test data. 
The second part of the script is the NNGP run. It will take some time. You can also start *run.R* on a cluster using *Rscript_NDVI.sh*. 
Now that the data is processed, you can also starts the runs for INLA with *INLA_stationary.R* and *Rscript_INLA.R*, for local GPs with *local_GPs.R* and *Rscript_local_GPs.sh*. 
You can run *inla_nonstat_noise.R* and *inla_nonstat_range.R*, but this will result in failure by RAM shortage - even if you use monstruous swap. 

## Compare the methods
Put the results in the */res* folder. 
Then run *aggregate_results.R*. 

## Plot the maps 
Run *maps.R*. The resulting maps will be in */maps_for_JCGS*

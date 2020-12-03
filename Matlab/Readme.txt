1) For running all of the following files, unzip All_electrodes_All_conds.zip into a folder and set this folder path in dataDir. 

2) Group_Elec_RASv2 and Electrode_Plotter_Group_McGurk_v3 are the files that you created with minor changes that I made to make it easier to run. You can change the analysis window in lines 27-30 to plot electrode activations accordingly. 

3) glm_rsquared: This essentially automates the process that we tried on excel by fitting a trend line and calculating the r-squared. This file fits a GLM and calculates the r-squared and standard errors and writes them out in an excel sheet which will then be used by analyse_rsquared and analyse_rsquared_indiv. Hence, glm_rsquared should be run before running analyse_squared and analyse_rsquared_indiv.

4) analyse_rsquared fits a glm to all the subjects together and calculates the r-squared and standard errors across all analysis windows and plots them out. 

5) analyse_rsquared_indiv plots the output of glm_squared in a lineplot-errorbar format. 


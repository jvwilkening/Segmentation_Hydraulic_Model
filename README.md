# Segmentation_Hydraulic_Model
Updated: July 2022

This repository contains the code developed and used in the manuscript "Exploring within-plant hydraulic trait variation: A test of the Vulnerability Segmentation Hypothesis".

This code is run in Python 3.10, and the necessary accompanying packages and versions are listed in the requirements.txt file. 

Author Info: Jean Wilkening (jvwilkening@berkeley.edu)

### Code
The code is run using the various "run..." scripts, one for each of the experiments and the case study as described in the manuscript. The parameter space for the experiments can be adjusted as denoted in the upper portion of the script files. The output from the scripts is then saved in an "Output" directory you will need to create. 

The plant hydraulics model is defined in SPAC_model.py and supporting functions for running the model are defined in model_functions.py. 

The file params_constants.py contains various constants that are used throughout the experiments, and the cleaned_XFT.txt, cleaned_XFT_Huber.txt, and cleaned_XFT_slope.txt contain data from the Xylem Functional Traits Database that is used for statistically generating parameter sets as described in the manuscript.

The "plot..." scripts contain the code to then retrieve the data generated from the run scripts, and plot the data as shown in manuscript. 

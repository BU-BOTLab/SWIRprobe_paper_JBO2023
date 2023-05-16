Instructions for running SWIR vs. NIR simulations from the paper titled "A shortwave infrared diffuse optical wearable probe for quantification of water and lipid content in emulsion phantoms using deep learning" (Spink et al., Journal of Biomedical Optics, 2023)

All scripts are run either in Matlab or Python. The Matlab version of Monte Carlo eXtreme (called MCXLAB) is required for the first script, titled createLUT.m. Furthermore, MCXLAB requires access to a GPU. Go to the following link for instructions on the download and installation of this software: http://mcx.space/wiki/index.cgi?Doc/MCXLAB. Step-by-step instructions for executing these simulations after Matlab, MCXLAB, and Python have all been installed are provided below.

1) Run createLUT.m in Matlab

The first step in the simulation process is to generate a look-up table (LUT) that maps diffuse reflectance to optical properties using Monte Carlo (MC) simulations. This LUT will be used in subsequent scripts for generating training and testing data for the deep neural networks (DNNs). Once the MC simulations are completed, check that the saved LUT is in same directory as the rest of the scripts for subsequent use (if, for whatever reason, it is not, then move the LUT to the same directory). You can modify the LUT name in the script, but by default it should save as "LUT_CW_multiDistance.mat".

2) Run generateTrainingData.m in Matlab 

This script will generate 100,000 sets of water and lipid concentrations (training data outputs), as well as their associated theoretical diffuse reflectance values for each of the 12 SWIR and NIR S-D pairs (training data inputs). This script will also call the "callLUT.m" function and the "chromophores_SWIR.txt" file. By default, this function should save the SWIR and NIR training data in the same directory in which the script is being run as two separate .mat files: R_SWIR_train_MC.mat and R_NIR_train_MC.mat. Each of these .mat files should contain the following variables:

R_SWIR(or R_NIR): 100,000 x 12 matrix of diffuse reflectance values
water: 100,000 x 1 array of water concentrations (as volume fractions)
lipid: 100,000 x 1 array of lipid concentrations (as volume fractions)
A_temp: 100,000 x 1 array of scattering amplitude values (in units of 1/mm)
b_temp: 100,000 x 1 array of scattering slope values (unitless)

3) Run generateTestingData.m in Matlab

This script will generate 25,000 set of water and lipid concentrations, as well as their associated theoretical diffuse reflectance values for each of the 12 SWIR and NIR S-D pairs (with noise added). Besides the size of the simulated dataset, the noise addition, and the names of the saved datasets, this step is identical to Step 2. The two .mat files will be named: R_SWIR_test_MC.mat and R_NIR_test_MC.mat.

4) Run swirDNN.py in Python

This script will construct and train a DNN with the parameters listed in the "SWIR vs. NIR comparison" section of the manuscript. This requires loading in the training data from Step 2, so be sure that the filepath and filenames called in the script match the saved file information from the previous step. The script will save the DNN as "DNN_SWIR.hdf5". It will also save the loss function history as "DNN_SWIR_loss.csv" and the normalization terms as "DNN_SWIR_norm_terms.csv".

WARNING: This repository also includes the previously constructed SWIR and NIR DNNs used in this paper in the folder "Preconstructed DNNs" with the same filenames. These will be overwritten if they are in the same directory as the swirDNN.py script.

5) Run nirDNN.py in Python

This step is identical to Step 4, but for the NIR DNN.

6) Run testSWIR.py in Python

This script will input the test SWIR test dataset to the SWIR DNN in order to estimate water and lipid concentrations. Error between the ground truth values and the DNN estimates will also be computed. Again, this script loads previously generated training data, testing data, and DNNs from Steps 2-4, so be sure that the saved file information matches the names and paths called in this script. The water and lipid estimates will be saved as CSV file called "DNN_SWIR_estimates.csv". 

NOTE: Training data and testing data will need to be generated using the steps described in this readme file, but the exact SWIR and NIR DNNs used in this paper are provided in the "Preconstructed DNNs" folder in this repository. If the user wants to load these DNNs directly instead of making/training them from scratch, Steps 4 and 5 can be skipped.

7) Run testNIR.py in Python
This step is identical to Step 6, but for the NIR DNN and test data.

8) Run simulationFigures.m in Matlab (optional)

This script will load the testing data, DNN estimates, and loss history files saved in previous steps and create the same visualizations shown in Figure 2 of the manuscript. 
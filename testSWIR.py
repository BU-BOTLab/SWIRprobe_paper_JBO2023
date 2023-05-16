
#%% Import 
import numpy as np # for performing mathematical operations
import scipy.io as spio # for loading .mat files
from tensorflow import keras # for implementing DNN
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Input
from tensorflow.keras.optimizers import Adam

#%% Load model
model = keras.models.load_model('./DNN_SWIR.hdf5') # load DNN 
matData_in = spio.loadmat('./R_SWIR_train_MC.mat') # load training data from mat file, for normalization of testing data
numInputs = 12 # 12 inputs refer to 12 S-D separation/wavelength combinations
numTestSamples = 25000 # set to length of testing data 
input_arr = matData_in['R_SWIR'] # 'R_SWIR' variable contains training data inputs

#%% Load test data and use model to estimate water and lipid concentrations
matData_test = spio.loadmat('./R_SWIR_test_MC.mat') # load testing data
eval_vec = matData_test['R_SWIR'] # 'R_SWIR' is matrix of test data inputs
eval_vec_norm = np.zeros((numTestSamples,numInputs)) # pre-allocate space for normalized test data inputs

for count in range(numInputs): # log normalization of test data
    eval_vec_norm[:,count] = (np.log10(eval_vec[:,count])-min(np.log10(input_arr[:,count])))/(max(np.log10(input_arr[:,count]))-min(np.log10(input_arr[:,count])))
 
p = model.predict(eval_vec_norm) # use DNN model to predict water and lipid estimates for test data inputs

groundtruth_e1=matData_test['water'] # define output variables for testing data , loaded from .mat file
groundtruth_e2=matData_test['lipid']
gt_e1=groundtruth_e1[:,0]
gt_e2=groundtruth_e2[:,0]
groundtruth_e = np.stack((gt_e1,gt_e2),axis=1) # concatenate output variables into one matrix

waterDiff = p[:,0] - groundtruth_e[:,0] # compute difference between DNN estimates and ground truth 
lipidDiff = p[:,1] - groundtruth_e[:,1] 

waterAE_mean = np.mean((100*waterDiff)) # compute mean and standard deviation of errors
lipidAE_mean = np.mean((100*lipidDiff))
waterAE_std = np.std((100*waterDiff))
lipidAE_std = np.std((100*lipidDiff))

# Display mean and standard deviation of errors (can also be done in subsequent Matlab script, described in readme file)
print(f"Absolute Error, Mean: {np.stack((waterAE_mean,lipidAE_mean),axis=0)}")
print(f"Absolute Error, St. Dev.: {np.stack((waterAE_std,lipidAE_std),axis=0)}")


#%%  Save estimates
np.savetxt('DNN_SWIR_estimates.csv', p, delimiter=',')

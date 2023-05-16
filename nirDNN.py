
#%% Import 
import numpy as np # for performing mathematical operations
import scipy.io as spio # for loading .mat files
from tensorflow import keras # for implementing DNN
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Input
from tensorflow.keras.optimizers import Adam

#%% Construct and train DNN
matData_in = spio.loadmat('./R_SWIR_train_MC.mat') # load training data from mat file
numInputs = 12 # 12 inputs refer to 12 S-D separation/wavelength combinations
numTrainSamples = 100000 # set to length of training data 
numHidden = 3 # set number of hidden layers in fully connected DNN
numNodes = 20 # set number of nodes/layer
numOutputs = 2 # set number of output variables (2 - water and lipid estimates)

input_arr = matData_in['R_NIR'] # 'R_NIR' variable contains training data inputs
input_arr_norm = np.zeros((numTrainSamples,numInputs)) # training data will need to be normalized; pre-allocate space for this
norm_terms = np.zeros((2,numInputs)) # the min and max values of training data will serve as normalization terms; pre-allocate space for these

for count in range(numInputs): # for each S-D pair, perform log normalization (min-max normalization of the log_10(inputs))
    input_arr_norm[:,count] = (np.log10(input_arr[:,count])-min(np.log10(input_arr[:,count])))/(max(np.log10(input_arr[:,count]))-min(np.log10(input_arr[:,count])))
    norm_terms[0,count] = max(np.log10(input_arr[:,count])) # store normalization terms 
    norm_terms[1,count] = min(np.log10(input_arr[:,count]))
    
groundtruth_1 = matData_in['water'] # define output variables for training data , loaded from .mat file
groundtruth_2 = matData_in['lipid']
gt_1 = groundtruth_1[:,0] 
gt_2 = groundtruth_2[:,0]
groundtruth = np.stack((gt_1,gt_2),axis=1) # concatenate ground truth water and lipid values into one matrix
    
model = Sequential() # define a linear neural network class
model.add(Input(shape=(numInputs,))) # add input layer to DNN
for l in range(numHidden): 
    model.add(Dense(numNodes,activation="relu")) # add hidden, fully-connected layers with ReLU activation functions
model.add(Dense(numOutputs, activation="linear")) # add output fully-connected layer with linear activation function
model.compile(loss="mse", optimizer=Adam(lr=.001)) # set loss function (mean squared error), training algorithm (Adam), and learning rate (0.001)

training = model.fit(x=input_arr_norm, y=groundtruth, epochs=3000, batch_size=256, validation_split=0.25) # train model; parse out 25% of training data to serve as validation set
loss_history = training.history["loss"] # store loss function history for future visualizations

#%% Save model, if desired
model.save('./DNN_NIR.hdf5') # save DNN
np.savetxt('DNN_NIR_norm_terms.csv', norm_terms, delimiter=',') # save normalization terms
np.savetxt('DNN_NIR_loss.csv', loss_history, delimiter=',') # save loss history
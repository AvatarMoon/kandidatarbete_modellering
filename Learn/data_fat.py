import scipy.io
import numpy as np
import os

# Skapar en väg till datand
data_paths = os.listdir('data_glucose_fat') 

for x in data_paths:
    
    data = scipy.io.loadmat('data_glucose_fat/' + x)
    data = data[x[0:-4]]

    concentration = data[0]['conc']

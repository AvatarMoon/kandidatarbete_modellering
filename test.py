import numpy as np 
import scipy.integrate as integrate 
from scipy.optimize import minimize
from scipy.optimize import Bounds
import matplotlib.pyplot as plt 
import math 
import pandas as pd
import os
import lhsmdu

def model1(t, x, k):
    k1, k2 = k
    
    dx1 = -k1*x[0] + k2*x[1] 
    dx2 = k1*x[0] - k2*x[1]
    
    return [dx1, dx2]

# Colour-blind friendly palette (use nice colors)
cb_palette1 = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

# get data 
data_G = pd.read_csv ("data_horses/glukos_FF_training.csv", sep=';')
data_I = pd.read_csv ("data_horses/insulin_FF_training.csv", sep=';')
data_G = data_G.sort_values(by=['time'])
data_I = data_I.sort_values(by=['time'])


# Extract times- and concentration-vectors
tG_vec = data_G['time'].values
tI_vec = data_I['time'].values  
cG_vec = data_G['conc'].values
cI_vec = data_I['conc'].values 

# Split the time-vectors at timepoint 20 minutes
index_G = np.searchsorted(tG_vec, 20)
index_I = np.searchsorted(tI_vec, 20)

tG_1 = tG_vec[0:index_G]
tG_2 = tG_vec[index_G:]

tI_1 = tI_vec[0:index_I]
tI_2 = tI_vec[index_I:]

x0 = [17,23]
b = [1,2]

first_sol_I = integrate.solve_ivp(model1, [0, 20], x0, method = "LSODA" , args = (b,), t_eval = tI_1)

x1 = first_sol_I.y[:,-1] + [0, 20]

print(x1)

second_sol_I = integrate.solve_ivp(model1, [20,240], x1 , method = "LSODA", args= (b, ), t_eval = tI_2)

sol_I = np.concatenate([first_sol_I.y, second_sol_I.y], axis = 1)

print(sol_I.shape)

plt.plot(tI_vec, sol_I[1])
plt.show()


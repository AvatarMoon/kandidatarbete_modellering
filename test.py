import numpy as np 
import scipy.integrate as integrate 
from scipy.optimize import minimize
from scipy.optimize import Bounds
import matplotlib.pyplot as plt 
import math 
import pandas as pd
import os
import lhsmdu

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

for n in tG_vec:
    if n <20
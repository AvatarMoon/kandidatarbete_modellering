import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import pandas as pd
 
data_G = pd.read_csv ("data_horses/glukos_FF_training.csv", sep=';')
data_I = pd.read_csv ("data_horses/insulin_FF_training.csv", sep=';')

# Extract times- and concentration-vectors
tG_vec = data_G['time'].values
tI_vec = data_I['time'].values  
cG_vec = data_G['conc'].values
cI_vec = data_I['conc'].values

plot1 = plt.figure(1)
plt.plot(tG_vec, cG_vec)
plt.title("Datapunkter glukos i plasma")
plt.show()

plot2 = plt.figure(2)
plt.plot(tI_vec, cI_vec)
plt.title("Datapunkter insulin i plasma")
plt.show()
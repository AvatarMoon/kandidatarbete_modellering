import matplotlib.pyplot as plt
import numpy as np 
import scipy.integrate as integrate 

# get data 
data_G = pd.read_csv ("data_horses/Glukos_new_FFaraber.csv", sep=';')
data_I = pd.read_csv ("data_horses/Insulin_new_FFaraber.csv", sep=';')
data_G = data_G.sort_values(by=['tid'])
data_I = data_G.sort_values(by=['tid'])


# Extract times- and concentration-vectors
tG_vec = data_G['tid'].values
tI_vec = data_I['tid'].values  

time_span = [tG_vec[0], tG_vec[-1]]

# Inital values for concentrations and optimal values for parameters 1 start guess
x0 = [5, 60, 34, 3, 2.5, 200]   
b = [1.65196917, 1.43044967, 2.19742316e+04, 5.59792273e+01, 2.84072466e+02, 6.46691412, 1.31869684e+01, 3.09662271e+03, 4.20835897e-02, 8.93446278e-02, 3.91588623, 1.04056980, 3.81839534, 2.73503274, 6.48068508, 1.15970185e+01, 3.11515367, 1.58857288e+01]

sol_qual = integrate.solve_ivp(open_loop, time_span, x0, method="LSODA", args=(b, ))

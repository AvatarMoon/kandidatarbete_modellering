import numpy as np 
import scipy.integrate as integrate 
from scipy.optimize import minimize
from scipy.optimize import Bounds
import matplotlib.pyplot as plt 
import math 
import pandas as pd
import os

# [Mörkblå, Gul, Ockraröd, Mörkgrön ,Olivgrön ,Ljusbeige (dålig), Ljusgult (dålig), Gul ]
cb_palette2 = ["#F4E3AF", "#F1CB6F", "#E16F3B", "#2D4E63", "#899964", "#F4E3AF", "#F1CB6F", "#E16F3B"]

# get data 
data_G = pd.read_csv ("data_horses/Glukos_new_FFaraber.csv", sep=';')
data_I = pd.read_csv ("data_horses/Insulin_new_FFaraber.csv", sep=';')
data_G = data_G.sort_values(by=['tid'])
data_I = data_I.sort_values(by=['tid'])


# Extract times- and concentration-vectors
tG_vec = data_G['tid'].values
tI_vec = data_I['tid'].values  
cG_vec = data_G['konc'].values
cI_vec = data_I['konc'].values 

Ge = 5 # konstant?

def open_loop(t,x,b): 
    c0, c1, c2, c3, b1, b2, b3, b4, b5, b10, b21, b22, b23, b25, b27, b100, f, v = b
     
    # Says that the concentrations can't be lower than zero 
    x[x < 0] = 0.0 

    # Concentrations in the model as input 
    G, I, E, C, M, H = x 

    # Glucose plasma [1]  G0 closed loop = 5mM 
    dG = f*b10*H/v + f*b5*C/v - b1*G - b3*I*G 

    # Insulin plasma [2]  60 pM 
    dI = b4*G - b2*I    

    # Glucacon plasma [3] E0 closed loop = 34 pM 
    dE = c0 + (c1/(c2 + I))*(Ge - G)*np.heaviside(Ge-G,1) - c3*E  

    # GLucose liver [4] C0 closed loop = 3 mmol 
    dC = b23 - b25*I - b22*G + b21*E - b5*C 

    # Glucose musle [5] M0 = 2.5 mmol 
    dM = 0.1*(v/f)*b3*G*I - b27*M 

    # Glucos intake [6]  H0 = 200 mmol 
    dH = - b100*H*G 

    return [dG, dI, dE, dC, dM, dH] 
 

def cost_function(b, yG_vec, yI_vec): 

    # Model parameters   
    x0 = [5, 60, 34, 3, 2.5, 200]  # initial closed loop   
    time_span_G = [tG_vec[0], tG_vec[-1]] 
    time_span_I = [tI_vec[0], tI_vec[-1]] 
    

    # Step 1: Solve ODE-system at points tG_vec
    sol_G = integrate.solve_ivp(open_loop, time_span_G, x0, method="LSODA", args=(b, ), t_eval=tG_vec) 
    sol_I = integrate.solve_ivp(open_loop, time_span_I, x0, method="LSODA", args=(b, ), t_eval=tI_vec) 

    
    # step 2: Solve ODE-system qualitative
    sol_qual = integrate.solve_ivp(open_loop, time_span_G, x0, method="LSODA", args=(b, ))

    G_model = sol_qual.y[0]
    I_model = sol_qual.y[1]
    E_model = sol_qual.y[2]
    C_model = sol_qual.y[3]
    M_model = sol_qual.y[4]

    # Step 3: Extract G and I model concentrations at t-points tG_vec and tI_vec
    yG_model = sol_G.y[0] 
    yI_model = sol_I.y[1] 

    # Step 4 : Build bounds for the concentrations and punnish the cost-func. if they go cross the bounds
    squared_sum = 0.0

    range_G = [4.5, 11] # mM 
    range_I = [38, 400] #pM 
    range_E = [28.68, 47.04] # pM 
    range_C = [0, 8] # mmol 
    range_M = [2, 13] # mmol 

    if any(G_model) > np.max(range_G):
       squared_sum += 100
    if any(G_model) < np.min(range_G):
        squared_sum += 100
    if any(I_model) > np.max(range_I):
       squared_sum += 100
    if any(I_model) < np.min(range_I):
        squared_sum += 100
    if any(E_model) > np.max(range_E):
       squared_sum += 100
    if any(E_model) < np.min(range_E):
        squared_sum += 100
    if any(C_model) > np.max(range_C):
       squared_sum += 100
    if any(C_model) < np.min(range_C):
        squared_sum += 100
    if any(M_model) > np.max(range_M):
       squared_sum += 100
    if any(M_model) < np.min(range_M):
        squared_sum += 100
    

    # Step 5: Calculate cost-function  
    squared_sum = np.sum((yG_model - yG_vec)**2+(yI_model -  yI_vec)**2) 

    return squared_sum 


# Note, a numerical optmizer require a starting guess for the parameters
start_values = [1.885, 198, 94, 0.0554, 0.0059, 0.1262, 0.00005, 0.4543, 0.185, 0.022, 0.00876, 0.0021, 0.08, 0.00026, 0.014, 0.3, 0.9, 15]
bound_low = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0.1])
bound_upp = np.repeat(np.inf, len(start_values))
bounds = Bounds(bound_low, bound_upp)
res = minimize(cost_function, start_values, method='Powell', args = (cG_vec, cI_vec), bounds=bounds) #lägg in constraints här 
 
# Print some statistics  
print("Optimal value found via Powells-method:") 
print(res.x) 
print("Value of cost-function") 
print(res.fun) 

# Plotting observed data at time-points 
data1 = plt.plot(tI_vec, cI_vec, color = cb_palette2[0]) 
data2 = plt.plot(tG_vec, cG_vec, color = cb_palette2[2])
# model1 =
# model2 = 
plt.legend(['Data insulin', 'Data glucose'], loc = 'upper left')
plt.title("Data") 

# Write the result to file
path_result_dir = "optimering/Bilder"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/model+data.pdf"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

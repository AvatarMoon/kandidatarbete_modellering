import numpy as np 
import scipy.integrate as integrate 
from scipy.optimize import minimize
from scipy.optimize import Bounds
import matplotlib.pyplot as plt 
import math 
import pandas as pd
import os
import lhsmdu

# [Mörkblå, Gul, Ockraröd, Mörkgrön ,Olivgrön ,Ljusbeige (dålig), Ljusgult (dålig), Gul ]
cb_palette2 = ["#F4E3AF", "#F1CB6F", "#E16F3B", "#2D4E63", "#899964", "#F4E3AF", "#F1CB6F", "#E16F3B"]

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

Ge = 5 # konstant?

def open_loop(t,x,b): 
    c0, c1, c2, c3, b1, b2, b3, b4, b5, b10, b21, b22, b23, b25, b27, b100, f, v = b
     
    # Says that the concentrations can't be lower than zero 
    x[x < 0] = 0.0 

    # Concentrations in the model as input 
    G, I, E, C, M, H = x 

    # Glucose plasma [1]  G0 closed loop = 5mM 
    dG = f*b10*H/(v+0.1) + f*b5*C/(v+0.1) - b1*G - b3*I*G 

    # Insulin plasma [2]  60 pM 
    dI = b4*G - b2*I    

    # Glucacon plasma [3] E0 closed loop = 34 pM 
    dE = c0 + (c1/(c2 + I))*(Ge - G)*np.heaviside(Ge-G,1) - c3*E  

    # GLucose liver [4] C0 closed loop = 3 mmol 
    dC = b23 - b25*I - b22*G + b21*E - b5*C 

    # Glucose musle [5] M0 = 2.5 mmol 
    dM = 0.1*(v/(f+0.1))*b3*G*I - b27*M 

    # Glucos intake [6]  H0 = 200 mmol 
    dH = - b100*H*G 

    return [dG, dI, dE, dC, dM, dH] 
 

def cost_function(b, yG_vec, yI_vec): 

    # Model parameters   
    x0 = [29, 67, 219, 19, 16, 1288]  # initial FF träning   
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

    range_G = [0, 50] # mM 
    range_I = [0, 5000] #pM 
    range_E = [0, 500] # pM 
    range_C = [0, 100] # mmol 
    range_M = [0, 140] # mmol 

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
    squared_sum = np.sum((yG_model - yG_vec))**2+np.sum((yI_model -  yI_vec)**2) 

    return squared_sum 

## Hypercube set up
randSeed = 2 # random number of choice
lhsmdu.setRandomSeed(randSeed) # Latin Hypercube Sampling with multi-dimensional uniformity
start = np.array(lhsmdu.sample(18, 4)) # Latin Hypercube Sampling with multi-dimensional uniformity (parameters, samples)

para, samples = start.shape

## intervals for the parameters
para_int = [0, 500]

minimum = (np.inf, None)

# Bounds for the model
bound_low = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0.1])
bound_upp = np.repeat(np.inf, para)
bounds = Bounds(bound_low, bound_upp)

start_values = [7.621, 800.505, 380.038, 2.24, 0.024, 0.510, 0.0002, 1.837, 0.748, 0.089, 0.035, 0.008, 0.323, 0.001, 0.057, 1.213, 3.64, 60.644]

for n in range(samples):
    c0 = start[0,n] * para_int[1]
    c1 = start[1,n] * para_int[1]
    c2 = start[2,n] * para_int[1]
    c3 = start[3,n] * para_int[1]
    b1 = start[4,n] * para_int[1]
    b2 = start[5,n] * para_int[1]
    b3 = start[6,n] * para_int[1]
    b4 = start[7,n] * para_int[1]
    b5 = start[8,n] * para_int[1]
    b10 = start[9,n] * para_int[1]
    b21 = start[10,n] * para_int[1]
    b22 = start[11,n] * para_int[1]
    b23 = start[12,n] * para_int[1]
    b25 = start[13,n] * para_int[1]
    b27 = start[14,n] * para_int[1]
    b100 = start[15,n] * para_int[1]
    f = start[16,n] * para_int[1]
    v = start[17,n] * para_int[1]
    
    res = minimize(cost_function, start_values, method='Powell', args = (cG_vec, cI_vec), bounds=bounds) #lägg in constraints här 

    if res.fun < minimum[0]:
        minimum = (res.fun, res.x)


# Print some statistics  
print("Optimal value found via Powells-method:") 
print(minimum[1]) 
print("Value of cost-function") 
print(minimum[0]) 

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


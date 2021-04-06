import numpy as np 
import scipy.integrate as integrate 
from scipy.optimize import minimize 
import matplotlib.pyplot as plt 
import math 
import pandas as pd
import os

# get data 
data_G = pd.read_csv ("data_horses/Glukos_new_FFaraber.csv", sep=';')
data_I = pd.read_csv ("data_horses/Insulin_new_FFaraber.csv", sep=';')
data_G = data_G.sort_values(by=['tid'])


global tG_vec 
global tI_vec 

# Extract times- and concentration-vectors
tG_vec = data_G['tid'].values
tI_vec = data_I['tid'].values  
cG_vec = data_G['konc'].values
cI_vec = data_I['konc'].values 


"""
# Tillfälliga parametrar, ev startgissningar, kan räknas om med ekv human vs horses
# Värdena är för människa från closed loop, ska dessa göras om till log? 
c0 = 1.8854 
c1 = 198 
c2 = 94 
c3 = 0.0554 
b1 = 0.0059 
b2 = 0.1262 
b3 = 0.00005 
b4 = 0.4543 
b5 = 0.185 
# b6= 0.8 tror inte de här är med i någon ODE
# b7= 0.8 
# b8= 0.8 
# b9 = 0.8 
b10 = 0.022 
b21 = 0.00876 
b22 = 0.0021 
b23 = 0.08 
b25 = 0.00026 
b27 = 0.014 
b100 = 0.3 
f = 0.9 
v = 15 
"""

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
    time_span_G = [tG_vec[0], tG_vec[-1]] #vad sätter vi här när det finns två olika?
    time_span_I = [tI_vec[0], tI_vec[-1]]

    # Step 1: Solve ODE-system  
    solG = integrate.solve_ivp(open_loop, time_span_G, x0, method="LSODA", args=(b, ), t_eval=tG_vec) 
    solI = integrate.solve_ivp(open_loop, time_span_I, x0, method="LSODA", args=(b, ), t_eval=tI_vec) 
    """class scipy.optimize.Bounds(lb, ub, keep_feasible=False) # vart ska denna in och hur sätter man olika bounds för olika konc. 
    """ 

    # Step 2: Extract G and I model concentrations 
    global yG_model
    global yI_model
    yG_model = solG.y[0] 
    yI_model = solI.y[1] 

    # for sats, om > värde, addera ngt till kostnadsfkn 
    """ range_G = [4.5, 11] # mM 
    range_I = [38, 400] #pM 
    range_E = [28.68, 47.04] # pM 
    range_C = [0, 8] # mmol 
    range_M = [2, 13] # mmol 
   

    for i in res: 
        if G > 11: 
            i += 100 
        elif I > 400: 
            i += 100 
        elif E > 47.04: 
            i+= 100 
        elif M > 13: 
            i+= 100 
        elif H > 200: 
            i+= 100  """   

    # Step 3: Calculate cost-function  
    squared_sum = np.sum((yG_model - cG_obs)**2+(yI_model -  cI_obs)**2) 

    return squared_sum 


# Note, a numerical optmizer require a starting guess, here I use the start-guess (0.022, 0.022, 0.022)
start_values = [1.885, 198, 94, 0.0554, 0.0059, 0.1262, 0.00005, 0.4543, 0.185, 0.022, 0.00876, 0.0021, 0.08, 0.00026, 0.014, 0.3, 0.9, 15]
res = minimize(cost_function, start_values, method='Powell', args = (cG_vec, cI_vec)) #lägg in constraints här 
"""loop för olika rates/startgissningar"""  
 

# Print some statistics  
print("Optimal value found via Powells-method:") 
print(res.x) 
print("Value of cost-function") 
print(res.fun) 

# Plotting observed data at time-points 0.1, ..., 2.0 (we have 50 data-points) 
plt.plot(t_vec, y_obs) 
plt.title("Simulated data and model") 
#plt.plot(t_vec, y_model) 
plt.show() 

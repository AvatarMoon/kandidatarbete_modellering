import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import math
import pandas as pd 

# get data
data_G = pd.read_csv ("C:/Users/bella/Documents/Studier/Kandidat/Data/Glukos_FFaraber.csv")
print("Glucose plasma data")
print(data_G)

data_I = pd.read_csv("C:/Users/bella/Documents/Studier/Kandidat/Data/Insulin_FFaraber.csv")

G_data = data_G.values # Turn the data into np-array (matrix)
I_data = data_I.values
y_obs = G_data[:, 1]
global t_data
global t_eval
tG_data = G_data[:, 0]
# tI_data = I_data[:, 0]
t_vec = G_data[1, :]

# Tillfälliga parametrar, ev startgissningar, kan räknas om med ekv human vs horses
c0 = 1.8854
c1 = 198
c2 = 94
c3 = 0.0554
b1 = 0.0059
b2 = 0.1262
b3 = 0.00005
b4 = 0.4543
b5 = 0.185
b6= 0.8
b7= 0.8
b8= 0.8
b9 = 0.8
b10 = 0.022
b21 = 0.00876
b22 = 0.0021
b23 = 0.08
b25 = 0.00026
b27 = 0.014
b100 = 0.3
Ge = 5
f = 0.9
v = 15

def open_loop(t,x,b):
    # b1, b3, b5, b6, b7, b8, b9, b10, f, v = b
    
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

def simulate_model(t_vec, sigma):
     
    # Model parameters 
    x0 = [5, 60, 34, 3, 2.5, 200] # closed loop initial
    # x0 = [1, 0.8, 0.5, 0.1, 0.2, 0.5] initial för hästar
    # t_vec borde hämtas från data
    time_span = [t_vec[0], t_vec[-1]] # 0 är första och -1 är sista
    rates = [b1, b3, b5, b6, b7, b8, b9, b10, f, v] #parametrar
    n_data_points = len(t_vec)
    
    # Note that by using t_eval I get the solution at the time-points in t_vec (I don't have to 
    # use any interpolation). 
    sol = integrate.solve_ivp(open_loop, time_span, x0, method="LSODA", args=(rates, ), t_eval=t_vec)
    
    
    # Extract x0 and add measurment noise 
    x0_obs = sol.y[0] + np.random.normal(loc=0.0, scale=sigma, size=n_data_points)
    # bara optimera en?
    return x0_obs
 
""" simulate data
t_vec = np.linspace(0, 240, num=20) #var 10e minut i ett dygn
# Set seed to reproduce
np.random.seed(123)
y0_obs = simulate_model(t_vec, 0.5) """

"""
    Cost-function for model1 which takes the model-parameters (k1, k2) and 
    observed data as input. 
"""

def cost_function(b, y_obs):
    
    # Model parameters  
    x0 = [5, 60, 34, 3, 2.5, 200]  # initial closed loop  
    time_span = [t_vec[0], t_vec[-1]]
    
    # Step 1: Solve ODE-system 
    sol = integrate.solve_ivp(open_loop, time_span, x0, method="LSODA", args=(b, ), t_eval=t_vec)
    """class scipy.optimize.Bounds(lb, ub, keep_feasible=False) # vart ska denna in och hur sätter man olika bounds för olika konc.
    """
    
    # Step 2: Extract x0 (simulated y-vec) , denna delen varierar beroende på del av data, def fkn här.
    #global y_model
    
    y0_model = sol.y[0]
    #y_model2 = sol.y[1]

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
    squared_sum = np.sum((y0_model - y_obs)**2)

    return squared_sum

# Note, a numerical optmizer require a starting guess, here I use the start-guess (0.022, 0.022, 0.022)
res = minimize(cost_function, [b1, b3, b5, b6, b7, b8, b9, b10, f, v], method='Powell', args = (y_obs, )) #lägg in constraints här
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
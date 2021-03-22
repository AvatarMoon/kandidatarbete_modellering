import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import math
import pandas as pd 

# get data
"""Glucose = pd.read_csv()
Insulin = pd.read_csv()"""

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
    # Glucose plasma [1]
    dG = f*b10*H/v + f*b5*C/v - b1*G - b3*I*G

    # Insulin plasma [2]
    dI = b4*G - b2*I   

    # Glucacon plasma [3]
    dE = c0 + (c1/(c2 + I))*(Ge - G)*np.heaviside(Ge-G,1) - c3*E 

    # GLucose liver [4]
    dC = b23 - b25*I - b22*G + b21*E - b5*C

    # Glucose musle [5]
    dM = 0.1*(v/f)*b3*G*I - b27*M

    # Glucos intake [6]
    dH = - b100*H*G

    return [dG, dI, dE, dC, dM, dH]

def simulate_model(t_vec, sigma):
     
    # Model parameters 
    x0 = [1, 0.8, 0.5, 0.1, 0.2, 0.5]   # initial
    t_vec = [0.231756251, 0.17399052, 3.205728606, 3.06, 6,229764595, 10,80577321, 15,76784945, 21,44814628, 27,10148577, 36,36133236, 44,19243988, 51.32843312, 59.14413644, 70.19568344, 81.86532375, 92.49999473, 96.041034, 98.8638527, 106.679556, 113.7770388, 120.9014789, 128.7595437, 135.6018945, 146.2904802, 157.7049885, 173.0080932, 194.1262815, 86.71186854, 84.17884126, 82.29471568, 80.7697004, 79.17921729, 79.10219631, 78.15291281]
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
 

# We simulate data at 144 data-points from 1 to 1440
t_vec = np.linspace(0, 240, num=20) #var 10e minut i ett dygn
# Set seed to reproduce
np.random.seed(123)
y0_obs = simulate_model(t_vec, 0.5)


"""
    Cost-function for model1 which takes the model-parameters (k1, k2) and 
    observed data as input. 
"""

def cost_function(b, y0_obs):
    
    # Model parameters  
    x0 = [1, 0.8, 0.5, 0.1, 0.2, 0.5]   # initial concentrations  
    time_span = [t_vec[0], t_vec[-1]]
    
    # Step 1: Solve ODE-system 
    sol = integrate.solve_ivp(open_loop, time_span, x0, method="LSODA", args=(b, ), t_eval=t_vec)
    """class scipy.optimize.Bounds(lb, ub, keep_feasible=False) # vart ska denna in och hur sätter man olika bounds för olika konc.
    """
    
    # Step 2: Extract x0 (simulated y-vec) , denna delen varierar beroende på del av data, def fkn här.
    #global y_model
    
    #y_model = sol.y[0] # byt ut siffran i y[], nr 0-13
    y0_model = sol.y[0]
    #y_model2 = sol.y[1]

    # for sats, om > värde, addera ngt till kostnadsfkn

    # Step 3: Calculate cost-function 
    squared_sum = np.sum((y0_model - y0_obs)**2)

    return squared_sum

# Note, a numerical optmizer require a starting guess, here I use the start-guess (0.022, 0.022, 0.022)
res = minimize(cost_function, [b1, b3, b5, b6, b7, b8, b9, b10, f, v], method='Powell', args = (y0_obs, )) #lägg in constraints här
 """loop för olika rates/startgissningar""" 

# Print some statistics 
print("Optimal value found via Powells-method:")
print(res.x)
print("Value of cost-function")
print(res.fun)

# Plotting observed data at time-points 0.1, ..., 2.0 (we have 50 data-points)
plt.plot(t_vec, y0_obs)
plt.title("Simulated data and model")
#plt.plot(t_vec, y_model)
plt.show()
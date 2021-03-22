import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import math
# hämta modell från cloosed loop
b1 = 0.0059
b3 = 0.00005
b5 = 0.185
b6 = 0.0102
b7 = 0.03
b8 = 0.022
b9 = 0.022
b10 = 0.022
f = 0.9
v = 15



def closed_loop(t,x):
    
    b1, b3, b5, b6, b7, b8, b9, b10, f, v = b
    x[x < 0] = 0.0
    S, L, G, C, I, W, E, M, H = x 

    # Moste of the regular parameters
   
    #l = 0.006
    #m = 0.04
    
    """ 
    k8 = 0.5275
    b13 = 0.0000095
    b14 = 0.0278
    # b17 = 0.02 # random siffra
    # b18 = 0.35
    # b19 = 0.004
    GLUT4 = 50.8472GtA = 135
    GLUT1 = 0.0283
    KmG4 = 146.851
    KmG1 = 1.082 
    Kgluc = 0.25
    Fat = 22
    #r = 0.04
    #p2u = 0.033
    # Ib =
    q1 = 0.0031
    q2 = 0.4054
    Gb = """

     

    # Stomach glucose [1]
    dS = b9*H-b8*S
    
    # Intestine glucose [2]
    dL = b8*S-b10*L
    
    # plasma glucose [3]  denna som vi ska matcha med data! Bör vi göra om MK fkn?
    dG = f*b10*L/v + f*b5*C/v - b1*G-b3*I*G
    
    #print("G = {}, W = {},I = {}".format(G, W, I))
    # plasma insulin [4]
    dI = b4*G + c*W*G-b2*I
    
    # plasma incretin [5]
    dW = b6*L - b7*W + s
    
    # plasma glucagon [6] # 
    # u = np.heaviside(Ge, G)  # lägg till u enligt artikel
    
    """if (Ge - G) >= 0:
        dE = c0 + (c1/(c2 + I*e))*(Ge - G) - c3*E 
    elif (Ge - G) < 0:
        dE = c0 - c3*E"""
    
    dE = c0 + (c1/(c2 + I*e))*(Ge - G) - c3*E             
    
    # liver glucose [7]
    dC = b23 - b25*I*e - b22*G + b21*E - b5*C
    
    # dynamics of glucose mass in muscle tissue [8]
    dM = 0.1*(v/f)*b3*G*I*e - b27*M

    dH = -b100*dS

    # Adipose tissue glucose mass (A) [9] Länk till cellulär nivå, ta bort om det går.
    #dA = k8*(GtA)/(KmG4 + GtA) + GLUT1*(GtA)/(KmG1 + GtA) - Kgluc*A

    # Dynamics of plasma leptin (Y) [10]
    #dY = b13*A*Fat - b14*Y
    #print("S = {}, Q = {}, m = {}, l = {}, I = {}".format(S, Q, m, l, I))
    
    # Dynamics of ghrelin concentration n plasma (Q) [11] RÄTTA DENNA!!!
    # dQ = b12*math.exp(-l*S)*math.exp(-m*I) - b11*Q

    # Linking the whole body model with the cellular one [13] RÄTTA DENNA!!!
    # dINSA = -INSA + (I-2)

    # Linking the whole body model with the cellular one [14]
    #dGtA = - q1*GtA + q2*(G - 2)

    return [dS, dL, dG, dI, dW, dE, dC, dM, dH]

    # Ranges 
    """ range_G = [4.5, 11] # mM
    range_I = [38, 400] #pM
    range_W = [5, 50]
    range_E = [28.68, 47.04]
    range_C = [0, 8]
    range_M = [2, 13]
    range_A = [30, 120]
    range_Y = [0, 0.6]
    range_Q = [8, 1146]
    # range S, L = none """
    
#startvärden
S0 = 4 #mmol
L0 = 14 #mmol
G0 = 5 #mmol
I0 = 60 #pM

t_vec = np.linspace(1, 220, num=20)
#time_span = [t_vec[0], t_vec[-1]] 

"""
    Function for simulating data for model1 at time-points provided by t-vec. The 
    data is simulated using an additivate Gaussian noise of strength sigma (sigma
    is provided by the user). Note, for model1 x2 is treated as the output. 
"""
def simulate_model(t_vec, sigma):
     
    # Model parameters 
    x0 = [S0, L0, G0, I0]   # initial
    time_span = [t_vec[0], t_vec[-1]] # 0 är första och -1 är sista
    rates = [b1, b3, b5, b6, b7, b8, b9, b10, f, v] #parametrar
    n_data_points = len(t_vec)
    
    # Note that by using t_eval I get the solution at the time-points in t_vec (I don't have to 
    # use any interpolation). 
    sol = integrate.solve_ivp(closed_loop, time_span, x0, method="LSODA", args=(rates, ), t_eval=t_vec)
    
    # Extract x0 and add measurment noise 
    x0_obs = sol.y[0] + np.random.normal(loc=0.0, scale=sigma, size=n_data_points)
    # bara optimera en?
    return x0_obs
 

# We simulate data at 144 data-points from 1 to 1440
t_vec = np.linspace(0.1, 1440, num=144) #var 10e minut i ett dygn
# Set seed to reproduce
np.random.seed(123)
y0_obs = simulate_model(t_vec, 0.5)


"""
    Cost-function for model1 which takes the model-parameters (k1, k2) and 
    observed data as input. 
"""

def cost_function(b, y0_obs):
    
    # Model parameters  
    x0 = [S0, L0, G0, I0]   # initial concentrations  
    time_span = [t_vec[0], t_vec[-1]]
    
    # Step 1: Solve ODE-system 
    sol = integrate.solve_ivp(closed_loop, time_span, x0, method="LSODA", args=(b, ), t_eval=t_vec)
    
    # Step 2: Extract x0 (simulated y-vec) , denna delen varierar beroende på del av data, def fkn här.
    global y_model
    
    #y_model = sol.y[0] # byt ut siffran i y[], nr 0-13
    y0_model = sol.y[0]
    #y_model2 = sol.y[1]

    # for sats, om > värde, addera ngt till kostnadsfkn

    # Step 3: Calculate cost-function 
    squared_sum = np.sum((y0_model - y0_obs)**2)

    return squared_sum

# Note, a numerical optmizer require a starting guess, here I use the start-guess (0.022, 0.022, 0.022)
res = minimize(cost_function, [b1, b3, b5, b6, b7, b8, b9, b10, f, v], method='Powell', args = (y0_obs, )) #lägg in constraints här

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

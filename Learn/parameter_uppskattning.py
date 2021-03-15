import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize
import matplotlib.pyplot as plt


"""                                                                                                                                                                      
     ODE-model of the reaction network:                                                                                                                                   
     x1 -> x2
     x2 -> x1                                                                                                                                                                          
    
     Args:                                                                                                                                                                
         x, vector of state-values (x1, x2)                                                                                                                           
         t, current time value                                                                                                                                            
         b, parameter-vector (must be passed as args, see below)                                                                                                          
    Returns:                                                                                                                                                             
         current value of the derivitative                                                                                                                                
"""
    # Startvärden
S0 = 4 #mmol
L0 = 14 #mmol
G0 = 5 #mmol
I0 = 60 #pM

def model1(t, x, b): 
    b1, b3, b5, b6, b7, b8, b9, b10, f, v = b
    H = 200 # tog H(0) så länge
    C = 20
    dx0 = b9*H-b8*x[0]  #S
    dx1 = b8*x[0]-b10*x[1] #L
    dx2 = f*b10*x[1]/v + f*b5*C/v - b1*x[2]-b3*x[3]*x[2] #G
    dx3 = b6*x[1] - b7*x[3] + x[0] # I
   
    return [dx0, dx1, dx2, dx3]

"""
    Function for simulating data for model1 at time-points provided by t-vec. The 
    data is simulated using an additivate Gaussian noise of strength sigma (sigma
    is provided by the user). Note, for model1 x2 is treated as the output. 
"""
def simulate_model1(t_vec, sigma):
     
    # Model parameters 
    x0 = [S0, L0, G0, I0]   # initial
    time_span = [t_vec[0], t_vec[-1]] # 0 är första och -1 är sista
    rates = [0.0059, 0.00005, 0.185, 0.0102, 0.03, 0.022, 0.022, 0.022, 0.9, 15] # vad innebär det här?
    n_data_points = len(t_vec)
    
    # Note that by using t_eval I get the solution at the time-points in t_vec (I don't have to 
    # use any interpolation). 
    sol = integrate.solve_ivp(model1, time_span, x0, method="LSODA", args=(rates, ), t_eval=t_vec)
    
    # Extract x0 and add measurment noise 
    x0_obs = sol.y[0] + np.random.normal(loc=0.0, scale=sigma, size=n_data_points)
    # bara optimera en?
    return x0_obs
 

# We simulate data at 144 data-points from 1 to 1440
t_vec = np.linspace(0.1, 1440, num=144) #var 10e minut i ett dygn
# Set seed to reproduce
np.random.seed(123)
y_obs = simulate_model1(t_vec, 0.5)


"""
    Cost-function for model1 which takes the model-parameters (k1, k2) and 
    observed data as input. 
"""

def cost_function(b, y_obs1):
    
    # Model parameters  
    x0 = [S0, L0, G0, I0]   # initial concentrations  
    time_span = [t_vec[0], t_vec[-1]]
    
    # Step 1: Solve ODE-system 
    sol = integrate.solve_ivp(model1, time_span, x0, method="LSODA", args=(b, ), t_eval=t_vec)
    
    # Step 2: Extract x0 (simulated y-vec) , denna delen varierar beroende på del av data, def fkn här.
    global y_model
    y_model = sol.y[0] # byt ut siffran i y[], nr 0-13

    # for sats, om > värde, addera ngt till kostnadsfkn

    # Step 3: Calculate cost-function 
    squared_sum = np.sum((y_model -  y_obs)**2)

    return squared_sum

# Note, a numerical optmizer require a starting guess, here I use the start-guess (0.022, 0.022, 0.022)
res = minimize(cost_function, [0.0059, 0.00005, 0.185, 0.0102, 0.03, 0.022, 0.022, 0.022, 0.9, 15], method='Powell', args = (y_obs, )) #lägg in constraints här

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

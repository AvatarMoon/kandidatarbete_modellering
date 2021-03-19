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
         k, parameter-vector (must be passed as args, see below)                                                                                                          
    Returns:                                                                                                                                                             
         current value of the derivitative                                                                                                                                
"""
def model(t, x, k):
    k1, k2 = k
    
    dx1 = -k1*x[0] + k2*x[1] 
    dx2 = k1*x[0] - k2*x[1]
    
    return [dx1, dx2]

"""
    Function for simulating data for model1 at time-points provided by t-vec. The 
    data is simulated using an additivate Gaussian noise of strength sigma (sigma
    is provided by the user). Note, for model1 x2 is treated as the output. 
"""
"""
def simulate_model1(t_vec, sigma):
     
    # Model parameters  
    x0 = [10.0, 20.0]   
    time_span = [t_vec[0], t_vec[-1]]
    rates = [1, 2]
    n_data_points = len(t_vec)
    
    # Note that by using t_eval I get the solution at the time-points in t_vec (I don't have to 
    # use any interpolation). 
    sol = integrate.solve_ivp(model, time_span, x0, method="LSODA", args=(rates, ), t_eval=t_vec)
    
    # Extract x2 and add measurment noise 
    x1_obs = sol.y[0] + np.random.normal(loc=0.0, scale=sigma, size=n_data_points)
    return x1_obs 
"""

def simulate_model2(t_vec, sigma): #vårt andra set av datapunkter
    # Model parameters  
    x0 = [10.0, 20.0]   
    time_span = [t_vec[0], t_vec[-1]]
    rates = [1, 2]
    n_data_points = len(t_vec)
    
    # Note that by using t_eval I get the solution at the time-points in t_vec (I don't have to 
    # use any interpolation). 
    sol = integrate.solve_ivp(model, time_span, x0, method="LSODA", args=(rates, ), t_eval=t_vec)
    
    # Extract x2 and add measurment noise 
    x2_obs = sol.y[1] + np.random.normal(loc=0.0, scale=sigma, size=n_data_points)
    return x2_obs


# We simulate data at 50 data-points from 0.1 to 2.0
#t_vec = np.linspace(0.1, 2.0, num=50
t_vec = [0.1    ,    0.11877551, 0.1255102 ,0.1332653 ,0.145510204, 0.287755,
 0.33965306 ,0.34142857 ,0.36020408 ,0.44897959 ,0.4877551  ,0.52653061,
 0.56530612 ,0.60408163 ,0.64285714, 0.68163265 ,0.72040816 ,0.75918367,
 0.79795918, 0.83673469 ,0.8755102 , 0.91428571 ,0.95306122 ,0.99183673,
 1.03061224 ,1.06938776 ,1.10816327, 1.14693878 ,1.18571429, 1.2244898,
 1.26326531, 1.30204082, 1.34081633, 1.37959184, 1.41836735, 1.45714286,
 1.49591837, 1.53469388, 1.57346939, 1.6122449 , 1.65102041, 1.68979592,
 1.72857143, 1.76734694 ,1.80612245, 1.84489796 ,1.88367347, 1.92244898,
 1.96122449, 2.        ]

# Set seed to reproduce
np.random.seed(123)
#y_obs1 = simulate_model1(t_vec, 0.5) #här sätter vi in riktiga datapunkter
y_obs1 = [ 9.4571847,  11.58919113, 12.21289911, 12.19099155, 13.42944331, 15.23699683,
 13.80917474, 15.35052827, 16.683327,   16.05016987, 16.52977931, 17.16621412,   
 18.26565248, 17.47295193, 17.81301694, 18.03340702, 19.54542822, 19.70676058,  
 19.26685423, 19.0924606,  19.38823203, 19.87314219, 18.7563838,  19.89791749,   
 18.75934244, 19.13566425, 19.96918621, 18.85538271, 19.54758238, 19.22925793,   
 19.57024161, 18.33232577, 18.87551843, 19.43765375, 20.27482079, 19.74526218,  
 19.85206147, 20.21120251, 19.4419978,  20.03679809, 19.50410396, 19.05339902,   
 19.73104368, 20.22163912, 20.11140778, 19.94272682, 21.15058481, 20.16595113,   
 20.45339015, 21.08711211]

y_obs2 = simulate_model2(t_vec, 0.5) 

print(y_obs2)

"""
    Cost-function for model1 which takes the model-parameters (k1, k2) and 
    observed data as input. 
"""

def cost_function(k, y_obs1, y_obs2):
    
    # Model parameters  
    x0 = [10.0, 20.0]   
    time_span = [t_vec[0], t_vec[-1]]
    
    # Step 1: Solve ODE-system 
    sol = integrate.solve_ivp(model, time_span, x0, method="LSODA", args=(k, ), t_eval=t_vec)
    
    # Step 2: Extract x2 (simulated y-vec)
    global y_model1  # global gör att variabeln finns utnför funktionen
    global y_model2 
    y_model1 = sol.y[0] 
    y_model2 = sol.y[1] 

    # Step 3: Calculate cost-function 
    squared_sum = np.sum((y_model1 -  y_obs1)**2+(y_model2 -  y_obs2)**2)

    return squared_sum

# Note, a numerical optmizer require a starting guess, here I use the start-guess (10.0, 10.0)
# (quite a bad start-guess)

res = minimize(cost_function, [10.0, 20.0], method='Powell', args = (y_obs1, y_obs2)) 

# Print some statistics 
print("Optimal value found via Powells-method:")
print(res.x)
print("True values")
print([1.0, 2.0])
print("Value of cost-function")
print(res.fun)


# Plotting observed data at time-points 0.1, ..., 2.0 (we have 50 data-points)
print("y = {}, t = {}".format(y_obs1, t_vec))
plot1 = plt.figure(1)
plt.plot(t_vec, y_obs1)
plt.title("Simulated data and model1")
plot1 = plt.figure(1)
plt.plot(t_vec, y_model1)

plot2 = plt.figure(2)
plt.plot(t_vec, y_obs2)
plt.title("Simulated data and model2")
plot2 = plt.figure(2)
plt.plot(t_vec, y_model2)
plt.show()

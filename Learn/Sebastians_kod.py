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
def model1(t, x, k):
    k1, k2 = k
    
    dx1 = -k1*x[0] + k2*x[1] 
    dx2 = k1*x[0] - k2*x[1]
    
    return [dx1, dx2]

"""
    Function for simulating data for model1 at time-points provided by t-vec. The 
    data is simulated using an additivate Gaussian noise of strength sigma (sigma
    is provided by the user). Note, for model1 x2 is treated as the output. 
"""
def simulate_model1(t_vec, sigma):
     
    # Model parameters  
    x0 = [10.0, 20.0]   
    time_span = [t_vec[0], t_vec[-1]]
    rates = [1, 2]
    n_data_points = len(t_vec)
    
    # Note that by using t_eval I get the solution at the time-points in t_vec (I don't have to 
    # use any interpolation). 
    sol = integrate.solve_ivp(model1, time_span, x0, method="LSODA", args=(rates, ), t_eval=t_vec)
    
    # Extract x2 and add measurment noise 
    x2_obs = sol.y[1] + np.random.normal(loc=0.0, scale=sigma, size=n_data_points)
    
    return x2_obs


# We simulate data at 50 data-points from 0.1 to 2.0
t_vec = np.linspace(0.1, 2.0, num=50)
# Set seed to reproduce
np.random.seed(123)
y_obs = simulate_model1(t_vec, 0.5)

# Plotting observed data at time-points 0.1, ..., 2.0 (we have 50 data-points)
plt.plot(t_vec, y_obs)
plt.title("Simulated data")
plt.show()

"""
    Cost-function for model1 which takes the model-parameters (k1, k2) and 
    observed data as input. 
"""

def cost_function(k, y_obs1):
    
    # Model parameters  
    x0 = [10.0, 20.0]   
    time_span = [t_vec[0], t_vec[-1]]
    
    # Step 1: Solve ODE-system 
    sol = integrate.solve_ivp(model1, time_span, x0, method="LSODA", args=(k, ), t_eval=t_vec)
    
    # Step 2: Extract x2 (simulated y-vec)
    y_model = sol.y[1] 

    # Step 3: Calculate cost-function 
    squared_sum = np.sum((y_model -  y_obs)**2)

    return squared_sum

# Note, a numerical optmizer require a starting guess, here I use the start-guess (10.0, 10.0)
# (quite a bad start-guess)
res = minimize(cost_function, [20.0, 20.0], method='Powell', args = (y_obs, ))

# Print some statistics 
print("Optimal value found via Powells-method:")
print(res.x)
print("True values")
print([1.0, 2.0])
print("Value of cost-function")
print(res.fun)
<<<<<<< Updated upstream
=======
plt.show()
>>>>>>> Stashed changes

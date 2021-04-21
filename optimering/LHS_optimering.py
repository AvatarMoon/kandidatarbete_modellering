import lhsmdu
import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize
import matplotlib.pyplot as plt

def model1(t, x, k):
    k1, k2 = k
    
    dx1 = -k1*x[0] + k2*x[1] 
    dx2 = k1*x[0] - k2*x[1]
    
    return [dx1, dx2]

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

def cost_function(k, y_obs):
    
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

randSeed = 2 # random number of choice
lhsmdu.setRandomSeed(randSeed) # Latin Hypercube Sampling with multi-dimensional uniformity
start = np.array(lhsmdu.sample(3, 20)) # Latin Hypercube Sampling with multi-dimensional uniformity

print(start)

res = minimize(cost_function, start, method='Powell', args = (y_obs, ))

print("Optimal value found via Powells-method:")
print(res.x)
print("True values")
print([1.0, 2.0])
print("Value of cost-function")
print(res.fun)
plt.show()


# fig = plt.figure()
# ax = fig.gca()
# ax.set_xticks(numpy.arange(0,1,0.1))
# ax.set_yticks(numpy.arange(0,1,0.1))

# plt.scatter(k[0], k[1], color='g', s = 30)
# plt.grid()
# plt.show()
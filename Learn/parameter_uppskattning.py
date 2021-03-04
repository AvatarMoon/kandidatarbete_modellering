import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize
import matplotlib.pyplot as plt

#dx1/dt = -k1*x1 + k2*x2
#dx2/dt = k1*x1 - k2*x2

def model(t, x, k):
    k1, k2 = k

    dx1 = -k1*x[0] + k2*x[1]
    dx2 = k1*x[0] - k2*x[1]

    return [dx1, dx2]

def simulate_model(t_vec, sigma):

    #model parameters 
    x0 = [10.0, 20.0]
    time_span = [t_vec[0], t_vec[-1]]
    rates = [1, 2]
    n_data_points = len(t_vec)

    sol = integrate.solve_ivp(model1, time_span, x0, method="LSODA", args = (rates, ), t_eval = t_vec)

    x2_obs = sol.y[1] + np.random.normal(loc=0.0, scale=sigma, size=n_data_points)

    return x2_obs

# We simulate data at 50 data-points from 0.1 to 2.0
t_vec = np.linspace(0.1, 0.2, num=50)
#set seed to reproduce
np.random.seed(123)
y_obs = simulate_model1(t_vec, 0.5)

#plotting observed data at time-points 0.1, ..., 2.0 (we have 50 data-points)
plt.plot(t_vec, y_obs)
plt.title("simulated data")

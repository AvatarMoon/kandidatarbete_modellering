import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt


"""
    ODE-model of the reaction network:
    phi1 -> x1 (phi1 = source term)
    x1 -> x2
    x2 -> x1
    x2 -> x3
    x3 -> phi2 (phi2 = breaks down)

    Adding rate constants the ODE is on the form 
    dx1/dt = k1 - k2*x1
    dx2/dt = k2*x1 - k3*x2
    dx3/dt = k3*x2 - k4*x3

    A typical ODE-model has the following args:
    Args:
        x, vector of state-values (x1, x2, x3)
        t, current time value 
        k, parameter-vector (must be passed as args, see below)
    Returns:
    	current value of the derivitative 
"""
def ode_model_example(t, x, k):

    # This is a neat way to access the elements of a vector 
    k1, k2, k3, k4 = k

    # Dynamics 
    dx1 = k1 - k2*x[0]
    dx2 = k2*x[0] - k3*x[1]
    dx3 = k3*x[1] - k4*x[2]
    
    return [dx1, dx2, dx3]


# Solving ODE-model 
time_span = [0.0, 4.0] # time-span for solving ODE 
rates = [1.0, 2.0, 3.0, 4.0] # parameter vector 
t_eval = np.linspace(0, 4, num=100)  # time-points to evaluate solution 
x0 = [1.0, 0.0, 0.0] # initial values 

# Solving ODE using LSODA-solver (a good solver), note rates enters via args and must be a tuple 
sol = integrate.solve_ivp(ode_model_example, time_span, x0, method="LSODA", args=(rates, ))

# Plot the solution 
plt.plot(sol.t, sol.y[0, :])
plt.plot(sol.t, sol.y[1, :])
plt.plot(sol.t, sol.y[2, :])
plt.legend(["x_1", "x_2", "x_3"])
plt.show()

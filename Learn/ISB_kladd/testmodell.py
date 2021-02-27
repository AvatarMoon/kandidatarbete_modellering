import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import pandas as pd

# Import experimental data
exemapledata = np.array(rawdata[1:],)
# Plot experimental data
# Define model 

def model(x1, x2):
    a = 3
    b = 15

    dG1 = a*x1 - b*x1**2 / x2
    dG2 = b*x2**2

    return [dG1, dG2]
# initial values [G1, G2]  
G0 = []  
# time
time = [0, 100]

sol = integrate.solve_ivp(model, time, G0, method='RK45')


# Initial guesses for parameters
# Range for parameters
# Parameter fitting




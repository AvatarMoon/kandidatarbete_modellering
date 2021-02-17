import numpy as np
import scipy.integrate as integrate
from scipy.integrate import RK45
import matplotlib.pyplot as plt

# Colour-blind friendly palette (use nice colors)
cb_palette = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

# function that returns dy/dt
def model(t,C):
    K_s = 1
    mu_max = .5
    Yield_S = 2
    mu = mu_max * (C[0]/(K_s+C[0]))

    dS = -Yield_S * mu * C[1]
    dX = mu * C[1]

    return [dS, dX]


# initial condition
C_start = [20, 0.5]

# time points
t = [0.0, 20.0]

# solve ODE
sol = integrate.solve_ivp(model,t,C_start, method = 'RK45')

# plot results
plt.plot(sol.t, sol.y[0, :], color = cb_palette[1])
plt.plot(sol.t, sol.y[1, :], color = cb_palette[0])
plt.xlabel('time')
plt.ylabel('Concentration')
plt.legend(["S", "X"])
plt.show()


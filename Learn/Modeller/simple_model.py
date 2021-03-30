import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

# Colour-blind friendly palette (use nice colors)
cb_palette = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

def model_normal(t,C):
    v = 15
    f = 0.9
    e = 1
    b1 = 0.0059     # /min
    b2 = 0.1262     # /min, insulin dissaperance rate
    b3 = 0.00005    # (pM*min)^-1
    b4 = 0.4543     # pM/mM * min
    b27 = .05       # /min
    G0 = 200        # mmol

    # Glucose in muscle
    dM = 0.1 * (v/f) * b3 * C[0] * C[1] * e - b27*C[2]

    # Insulin in plasma
    dI = b4 * C[0] - b2*C[1]

    # Glucose in plasma
    dG = G0 - b1 * C[0] - b3 * C[0] * C[1] 

    return [dG, dI, dM]

# inital values fasting [G, I, M]
C0 = [5e-3, 60e-12, 2.5e3]

# time
time = [0, 100]

sol = integrate.solve_ivp(model_normal,time,C0, method= 'LSODA')

#plottar

plt.plot(sol.t, sol.y[0, :], color = cb_palette[0])
plt.plot(sol.t, sol.y[1, :], color = cb_palette[1])
plt.plot(sol.t, sol.y[2, :], color = cb_palette[2])
plt.xlabel('time')
plt.ylabel('Concentration')
plt.legend(["G", "I", "M"])
plt.show()
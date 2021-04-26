import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import math

# Colour-blind friendly palette (use nice colors)
cb_palette = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

# [Mörkblå, Gul, Ockraröd, Mörkgrön ,Olivgrön ,Ljusbeige (dålig), Ljusgult (dålig), Gul ]
cb_palette2 = ["#F4E3AF", "#F1CB6F", "#E16F3B", "#2D4E63", "#899964", "#F4E3AF", "#F1CB6F", "#E16F3B"]

def open_loop(t,x):
    
    # Says that the concentrations can't be lower than zero
    x[x < 0] = 0.0

    # Concentrations in the model as input
    G, I, C, M, H = x

    # Parametrar
    k1 = 0.1
    k2 = 0.1
    k3 = 0.1
    k4 = 0.1
    k5 = 0.1 

    
    L= 500 # Startvärde glukos i levern

    # Glucose plasma [1]
    dG = k4*I*C + k1*H - k2*M

    # Insulin plasma [2]
    dI = -k4*I*C + k3*G

    # GLucose liver [4]
    dC = -k4*I*C + L

    # Glucose musle [5]
    dM = k2*G - k5*M

    # Glucose intake 
    dH = -k1 * H

    return [dG, dI, dC, dM, dH]

# Time span
t_vec = np.linspace(0, 240, num=20)
time_span = [t_vec[0], t_vec[-1]]

# Initial values to model
x0 = [60, 5, 34, 3, 70]  # G, I, C, M

# Solve ODE
sol = integrate.solve_ivp(open_loop, time_span, x0, method="LSODA", t_eval=t_vec)

# Plot model

ymodel = sol.y[0]
plt.plot(t_vec, ymodel)
plt.title('Open loop-model')
plt.show()
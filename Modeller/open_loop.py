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
    G, I, E, C, M, H = x

    #c = 0.1060
    c0 = 1.8854
    c1 = 198
    c2 = 94
    c3 = 0.0554
    b1 = 0.0059
    b2 = 0.1262
    b3 = 0.00005
    b4 = 0.4543
    b5 = 0.185
    b10 = 0.022
    b21 = 0.00876
    b22 = 0.0021
    b23 = 0.08
    b25 = 0.00026
    b27 = 0.014
    b100 = 0.3
    Ge = 5
    f = 0.9
    v = 15


    # Glucose plasma [1]
    dG = f*b10*H/v + f*b5*C/v - b1*G - b3*I*G

    # Insulin plasma [2]
    dI = b4*G - b2*I   

    # Glucacon plasma [3]
    dE = c0 + (c1/(c2 + I))*(Ge - G)*np.heaviside(Ge-G,1) - c3*E 

    # GLucose liver [4]
    dC = b23 - b25*I - b22*G + b21*E - b5*C

    # Glucose musle [5]
    dM = 0.1*(v/f)*b3*G*I - b27*M

    # Glucos intake [6]
    dH = -b100*G*H

    return [dG, dI, dE, dC, dM, dH]

# Time span
t_vec = np.linspace(0, 240, num=20)
time_span = [t_vec[0], t_vec[-1]]

# Initial values to model
x0 = [5, 60, 34, 3, 2.5, 200]

# Solve ODE
sol = integrate.solve_ivp(open_loop, time_span, x0, method="LSODA", t_eval=t_vec)

# Plot model

ymodel = sol.y[1]
plt.plot(t_vec, ymodel)
plt.title('Open loop-model')
plt.show()
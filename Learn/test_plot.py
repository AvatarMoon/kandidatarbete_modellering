import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import scipy.integrate as integrate

# Colour-blind friendly palette (use nice colors)
cb_palette1 = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

# [Mörkblå, Gul, Ockraröd, Mörkgrön ,Olivgrön ,Ljusbeige (dålig), Ljusgult (dålig), Gul ]
cb_palette2 = ["#F4E3AF", "#F1CB6F", "#E16F3B", "#2D4E63", "#899964", "#F4E3AF", "#F1CB6F", "#E16F3B"]

# Hämta data
data_G = pd.read_csv ('data_horses/Glukos_new_FFaraber.csv', sep=';')
data_I = pd.read_csv ('data_horses/insulin_new_FFaraber.csv', sep=';')
data_G = data_G.sort_values(by = ['tid'])
data_I = data_I.sort_values(by = ['tid'])

# Tidsvektorer
tG_vec = data_G['tid'].values
tI_vec = data_I['tid'].values
cG_vec = data_G['konc'].values
cI_vec = data_I['conc'].values


time_span_G = [tG_vec[0], tG_vec[-1]]
time_span_I = [tI_vec[0], tI_vec[-1]]

# Modellen
def open_loop(t,x,b): 
    c0, c1, c2, c3, b1, b2, b3, b4, b5, b10, b21, b22, b23, b25, b27, b100, f, v = b

    Ge = 5
     
    # Says that the concentrations can't be lower than zero 
    x[x < 0] = 0.0 

    # Concentrations in the model as input 
    G, I, E, C, M, H = x 

    # Glucose plasma [1]  G0 closed loop = 5mM 
    dG = f*b10*H/v + f*b5*C/v - b1*G - b3*I*G 

    # Insulin plasma [2]  60 pM 
    dI = b4*G - b2*I    

    # Glucacon plasma [3] E0 closed loop = 34 pM 
    dE = c0 + (c1/(c2 + I))*(Ge - G)*np.heaviside(Ge-G,1) - c3*E  

    # GLucose liver [4] C0 closed loop = 3 mmol 
    dC = b23 - b25*I - b22*G + b21*E - b5*C 

    # Glucose musle [5] M0 = 2.5 mmol 
    dM = 0.1*(v/f)*b3*G*I - b27*M 

    # Glucos intake [6]  H0 = 200 mmol 
    dH = - b100*H*G 

    return [dG, dI, dE, dC, dM, dH] 

# Inital values for concentrations and optimal values for parameters 1 start guess
x0 = [5, 60, 34, 3, 2.5, 200]   
b = [1.65196917, 1.43044967, 2.19742316e+04, 5.59792273e+01, 2.84072466e+02, 6.46691412, 1.31869684e+01, 3.09662271e+03, 4.20835897e-02, 8.93446278e-02, 3.91588623, 1.04056980, 3.81839534, 2.73503274, 6.48068508, 1.15970185e+01, 3.11515367, 1.58857288e+01]

sol_qual_G = integrate.solve_ivp(open_loop, time_span_G, x0, method="LSODA", args=(b, ), t_eval = tG_vec)
sol_qual_I = integrate.solve_ivp(open_loop, time_span_I, x0, method="LSODA", args=(b, ), t_eval = tI_vec)

G_model = sol_qual_G.y[0]
I_model = sol_qual_I.y[1]
E_model = sol_qual_G.y[2]
C_model = sol_qual_G.y[3]
M_model = sol_qual_G.y[4]

# Constrains glukos (G)
xG_coordinates = [tG_vec[0],tG_vec[-1]]
yG1_coordinates = [4.5,4.5] #mM (human)
yG2_coordinates = [11,11]   #mM (human)

# Constrains insulin (I)
xI_coordinates = [tI_vec[0],tI_vec[-1]]
yI1_coordinates = [38,38]   #pM (human)
yI2_coordinates = [400,400] #pM (human)

# Constrains glucagon (E)
xE_coordinates = [tG_vec[0],tG_vec[-1]] 
yE1_coordinates = [28.68,28.68]  #pM (human)
yE2_coordinates = [47.04,47.04]  #pM (human)
 
 # Constrains Glukos i lever (C)
xC_coordinates = [tG_vec[0],tG_vec[-1]] 
yC1_coordinates = [0,0]  # mmol (human)
yC2_coordinates = [8,8]  # mmol (human)

 # Constrains Glukos i muskel (M)
xM_coordinates = [tG_vec[0],tG_vec[-1]] 
yM1_coordinates = [2,2]    # mmol (human)
yM2_coordinates = [13,13]  # mmol (human)

# plotta glukos
lw = 2.0
plot1 = plt.figure(1)
frame1 = plot1.add_axes((.1, .3, 8., .6))
line1, = plt.plot(xG_coordinates, yG1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xG_coordinates, yG2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(data_G['tid'].values, data_G['konc'].values, label = 'Glukos', linestyle="-", linewidth=lw, color=cb_palette1[2])
line4, = plt.plot(data_G['tid'].values, G_model, label = 'Glukos', linestyle="-", linewidth=lw, color=cb_palette1[4])
plt.legend((line1, line2, line3, line4), ("Lägsta gräns", "Högsta gräns", "Data","modell"))
difference = cG_vec - G_model
frame2 = plot1.add_axes((.1, .1, .8, .2))
plt.plot(tG_vec, difference, 'or')
plt.xlabel("time", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)

plt.title("Glucose in plasma")

plt.show()

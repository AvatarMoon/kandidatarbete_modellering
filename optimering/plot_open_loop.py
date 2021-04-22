import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import scipy.integrate as integrate

# Colour-blind friendly palette (use nice colors)
cb_palette1 = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

# Hämta data
data_G = pd.read_csv ('data_horses\glukos_FF_träning.csv', sep=';')
data_I = pd.read_csv ('data_horses\insulin_FF_träning.csv', sep=';')
data_G = data_G.sort_values(by = ['time'])
data_I = data_I.sort_values(by = ['time'])

# Tidsvektorer
tG_vec = data_G['time'].values
tI_vec = data_I['time'].values

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
x0 = [29, 67, 219, 19, 16, 1288]   
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
yG1_coordinates = [0,0] #mM (human)
yG2_coordinates = [50,50]   #mM (human)

# Constrains insulin (I)
xI_coordinates = [tI_vec[0],tI_vec[-1]]
yI1_coordinates = [0,0]   #pM (human)
yI2_coordinates = [5000,5000] #pM (human)

# Constrains glucagon (E)
xE_coordinates = [tG_vec[0],tG_vec[-1]] 
yE1_coordinates = [0,0]  #pM (human)
yE2_coordinates = [500,500]  #pM (human)
 
 # Constrains glukos i lever (C)
xC_coordinates = [tG_vec[0],tG_vec[-1]] 
yC1_coordinates = [0,0]  # mmol (human)
yC2_coordinates = [100,100]  # mmol (human)

 # Constrains glukos i muskel (M)
xM_coordinates = [tG_vec[0],tG_vec[-1]] 
yM1_coordinates = [0,0]    # mmol (human)
yM2_coordinates = [140,140]  # mmol (human)

# plotta glukos
lw = 2.0
plot1 = plt.figure(1)
line1, = plt.plot(xG_coordinates, yG1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xG_coordinates, yG2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(data_G['time'].values, data_G['conc'].values, label = 'Glukos', linestyle="-", linewidth=lw, color=cb_palette1[7])
line4, = plt.plot(data_G['time'].values, G_model, label = 'Glukos', linestyle="-", linewidth=lw, color=cb_palette1[5])
plt.legend((line4, line3, line2, line1), ("Modell", "Data", "Högsta gräns","Lägsta gräns"))
plt.xlabel("time", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Glucose in plasma")

# Sparar figur i plot constrains, glukos
# Write the result to file
path_result_dir = "optimering/Bilder/plot_constrains"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_glucose.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta insulin
lw = 2.0
plot1 = plt.figure(2)
line1, = plt.plot(xI_coordinates, yI1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xI_coordinates, yI2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(data_I['time'].values, data_I['conc'].values, label = 'Insulin', linestyle="-", linewidth=lw, color=cb_palette1[7])
line4, = plt.plot(data_I['time'].values, I_model, label = 'Insulin', linestyle="-", linewidth=lw, color=cb_palette1[5])
plt.legend((line4, line3, line2, line1), ("Modell", "Data", "Högsta gräns","Lägsta gräns"))
plt.xlabel("time", fontsize=12), plt.ylabel("Insulin koncentration", fontsize=12)
plt.title("Insulin in plasma")

# Sparar figur i plot constrains, insulin
# Write the result to file
path_result_dir = "optimering/Bilder/plot_constrains"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_insulin.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta glukagon
lw = 2.0
plot1 = plt.figure(3)
line1, = plt.plot(xE_coordinates, yE1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xE_coordinates, yE2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(data_G['time'].values, E_model, label = 'Glucagon', linestyle="-", linewidth=lw, color=cb_palette1[5]) # Lägga till modellen
plt.legend((line3, line2, line1), ("Modell", "Högsta gräns", "Lägsta gräns"))
plt.xlabel("time", fontsize=12), plt.ylabel("Glukagon koncentration", fontsize=12)
plt.title("Glukagon i plasman")

# Sparar figur i plot constrains, glukagon
# Write the result to file
path_result_dir = "optimering/Bilder/plot_constrains"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_glukagon.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta Glukos i lever 
lw = 2.0
plot1 = plt.figure(4)
line1, = plt.plot(xC_coordinates, yC1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xC_coordinates, yC2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(data_G['time'].values, C_model, label = 'Glukos i levern', linestyle="-", linewidth=lw, color=cb_palette1[5]) # Lägga till modellen
plt.legend((line3, line2, line1), ("Modell", "Högsta gräns", "Lägsta gräns"))
plt.xlabel("time", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Glukos i levern")

# Sparar figur i plot constrains, glukos i levern
# Write the result to file
path_result_dir = "optimering/Bilder/plot_constrains"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_glukoslevern.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta Glukos i muskeln 
lw = 2.0
plot1 = plt.figure(5)
line1, = plt.plot(xM_coordinates, yM1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xM_coordinates, yM2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(data_G['time'].values, M_model, label = 'Glukos i muskeln', linestyle="-", linewidth=lw, color=cb_palette1[5]) # Lägga till modellen
plt.legend((line3, line2, line1), ("Modell", "Högsta gräns", "Lägsta gräns"))
plt.xlabel("time", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Glukos i muskeln")
plt.show()

# Sparar figur i plot constrains, glukos i muskeln
# Write the result to file
path_result_dir = "optimering/Bilder/plot_constrains"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_glukosmuskeln.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

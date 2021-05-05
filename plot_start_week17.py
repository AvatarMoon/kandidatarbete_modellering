import numpy as np 
import scipy.integrate as integrate 
from scipy.optimize import minimize
from scipy.optimize import Bounds
import matplotlib.pyplot as plt 
import math 
import pandas as pd
import os
import lhsmdu
from tqdm import tqdm
import copy
import datetime

# Colour-blind friendly palette (use nice colors)
cb_palette1 = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

# get data 
data_G = pd.read_csv ("data_horses/glukos_FF_training.csv", sep=';')
data_I = pd.read_csv ("data_horses/insulin_FF_training.csv", sep=';')
data_G = data_G.sort_values(by=['time'])
data_I = data_I.sort_values(by=['time'])


# Extract times- and concentration-vectors
tG_vec = data_G['time'].values
tI_vec = data_I['time'].values  
cG_vec = data_G['conc'].values
cI_vec = data_I['conc'].values 

# Split the time-vectors at timepoint 20 minutes
index_G = np.searchsorted(tG_vec, 20)
index_I = np.searchsorted(tI_vec, 20)

tG_1 = tG_vec[0:index_G]
tG_2 = tG_vec[index_G:]

tI_1 = tI_vec[0:index_I]
tI_2 = tI_vec[index_I:]

time_span1 = [0, 20]

time_span_G2 = [20, tG_2[-1]]
time_span_I2 = [20, tI_2[-1]]


def open_loop(t,x,b): 
    # Parameters for the models as input
    k1, k2, k3, k4, k5, k6, k7, k8 = b
     
    # # Says that the concentrations can't be lower than zero 
    # x = copy.deepcopy(x_old) # can't modify the input so we copy it

    # x[x < 0] = 0.0  

    # Scaling factor for units
    scal_factor = 1e-9

    # Concentrations in the model as input 
    G, I, C, M, H, E, F = x

    # Glucose plasma [1]
    dG = k4*C + k1*H - k2*G*I

    # Insulin plasma [2]
    dI = k3*G - k2*I*G*scal_factor

    # GLucose liver [3]
    dC = -k4*C + k6*scal_factor*E + k7*F

    # Glucose musle [4]
    dM = k2*scal_factor*G*I - k5*M

    # Glucose intake [5]
    dH = -k1*H

    # Glucagon in plasma [6]
    dE = k8 - k2*E*G 

    # Fettreserve [7]
    dF = -k7*F

    return [dG, dI, dC, dM, dH, dE, dF]

b = [0.0119979 , 0.00333303, 0.00152618, 0.00234708, 0.00571778, 1.07024714, 0.00260788, 0.27374058]
x0 = [30, 100, 34, 60, 70, 50, 400]  # G, I, C, M, H, E, F 
time_span_G = [tG_vec[0], tG_vec[-1]] 
time_span_I = [tI_vec[0], tI_vec[-1]] 

#Injection
inj = 2742

# Solve ODE-system qualitative
first_sol_qual = integrate.solve_ivp(open_loop, [0,20], x0, method="Radau", args=(b, ))

# Simulate the injection
x2 = first_sol_qual.y[:, -1] + [0, inj, 0, 0, 0, 0, 0]

# Solve ODE-system after 20 miunutes with injection
second_sol_qual = integrate.solve_ivp(open_loop, [20,240], x2, method = "Radau", args = (b, ))

sol_qual = np.concatenate([first_sol_qual.y, second_sol_qual.y], axis = 1)
time_span = np.concatenate([first_sol_qual.t, second_sol_qual.t])

G_model = sol_qual[0]
I_model = sol_qual[1]
C_model = sol_qual[2]
M_model = sol_qual[3]
H_model = sol_qual[4]
E_model = sol_qual[5]
F_model = sol_qual[6]

# Time span
xT_coordinates = [tG_vec[0],tG_vec[-1]]

# Constrains glukos (G)
yG1_coordinates = [0,0] #mM (human)
yG2_coordinates = [50,50]   #mM (human)

# Constrains insulin (I)
yI1_coordinates = [0,0]   #pM (human)
yI2_coordinates = [5000,5000] #pM (human)
 
 # Constrains glukos i lever (C) 
yC1_coordinates = [0,0]  # mmol (human)
yC2_coordinates = [100,100]  # mmol (human)

 # Constrains glukos i muskel (M) 
yM1_coordinates = [0,0]    # mmol (human)
yM2_coordinates = [140,140]  # mmol (human)

 # Constrains glucose intake (H)
yH1_coordinates = [0,0]    
yH2_coordinates = [500,500]  

 # Constrains glucagon in plasma (E)
yE1_coordinates = [0,0]    
yE2_coordinates = [500,500] 

 # Constrains glucagon in plasma (E)
yF1_coordinates = [0,0]    
yF2_coordinates = [500,500] 


# plotta glukos
lw = 2.0
plot1 = plt.figure(1)
# G_plot = plt.subplot(121)
line1, = plt.plot(xT_coordinates, yG1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xT_coordinates, yG2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(data_G['time'].values, data_G['conc'].values, label = 'Glukos', linestyle="-", linewidth=lw, color=cb_palette1[7])
line4, = plt.plot(time_span, G_model, label = 'Glukos', linestyle="-", linewidth=lw, color=cb_palette1[5])
plt.legend((line4, line3, line2, line1), ("Modell", "Data", "Högsta gräns","Lägsta gräns"))
plt.xlabel("time", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Glucose in plasma")

# # Residual plot for glucose
# G_res = plt.subplot(122)
# difference_G = cG_vec - G_model
# plt.scatter(difference_G, G_model, s = 10 , color = cb_palette1[1])

# Sparar figur i plot constrains, glukos
# Write the result to file
path_result_dir = "optimering/Bilder/test_plot_week17"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.makedirs(path_result_dir, exist_ok=True)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_glucose.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta insulin
lw = 2.0
plot1 = plt.figure(2)
# I_plot = plt.subfig(121)
line1, = plt.plot(xT_coordinates, yI1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xT_coordinates, yI2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(data_I['time'].values, data_I['conc'].values, label = 'Insulin', linestyle="-", linewidth=lw, color=cb_palette1[7])
line4, = plt.plot(time_span, I_model, label = 'Insulin', linestyle="-", linewidth=lw, color=cb_palette1[5])
plt.legend((line4, line3, line2, line1), ("Modell", "Data", "Högsta gräns","Lägsta gräns"))
plt.xlabel("time", fontsize=12), plt.ylabel("Insulin koncentration", fontsize=12)
plt.title("Insulin i plasman")

# # Residual plot for insulin
# I_res = plt.subplot(122)
# difference_I = cI_vec - I_model
# plt.scatter(difference_I, I_model, s = 10 , color = cb_palette1[1])

# Sparar figur i plot constrains, insulin
# Write the result to file
path_result_dir = "optimering/Bilder/test_plot_week17"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_insulin.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta Glukos i lever 
lw = 2.0
plot1 = plt.figure(3)
line1, = plt.plot(xT_coordinates, yC1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xT_coordinates, yC2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(time_span, C_model, label = 'Glukos i levern', linestyle="-", linewidth=lw, color=cb_palette1[5]) # Lägga till modellen
plt.legend((line3, line2, line1), ("Modell", "Högsta gräns", "Lägsta gräns"))
plt.xlabel("time", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Glukos i levern")

# Sparar figur i plot constrains, glukos i levern
# Write the result to file
path_result_dir = "optimering/Bilder/test_plot_week17"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_glukoslevern.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta Glukos i muskeln 
lw = 2.0
plot1 = plt.figure(4)
line1, = plt.plot(xT_coordinates, yM1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xT_coordinates, yM2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(time_span, M_model, label = 'Glukos i muskeln', linestyle="-", linewidth=lw, color=cb_palette1[5]) # Lägga till modellen
plt.legend((line3, line2, line1), ("Modell", "Högsta gräns", "Lägsta gräns"))
plt.xlabel("time", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Glukos i muskeln")


# Sparar figur i plot constrains, glukos i muskeln
# Write the result to file
path_result_dir = "optimering/Bilder/test_plot_week17"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_glukosmuskeln.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)


# plotta Glukos inktake 
lw = 2.0
plot1 = plt.figure(5)
line1, = plt.plot(xT_coordinates, yH1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xT_coordinates, yH2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(time_span, H_model, label = 'Glukos intag', linestyle="-", linewidth=lw, color=cb_palette1[5]) # Lägga till modellen
plt.legend((line3, line2, line1), ("Modell", "Högsta gräns", "Lägsta gräns"))
plt.xlabel("time", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Glukos intag")


# Sparar figur i plot constrains, glukos i muskeln
# Write the result to file
path_result_dir = "optimering/Bilder/test_plot_week17"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_glukos_intag.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta Glucagon in plasma
lw = 2.0
plot1 = plt.figure(6)
line1, = plt.plot(xT_coordinates, yE1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xT_coordinates, yE2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(time_span, E_model, label = 'Glukagon i plasma', linestyle="-", linewidth=lw, color=cb_palette1[5]) # Lägga till modellen
plt.legend((line3, line2, line1), ("Modell", "Högsta gräns", "Lägsta gräns"))
plt.xlabel("time", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Glukagon i plasma")


# Sparar figur i plot constrains, glukos i muskeln
# Write the result to file
path_result_dir = "optimering/Bilder/test_plot_week17"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_glucagon_plasma.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta Fettreserve 
lw = 2.0
plot1 = plt.figure(7)
line1, = plt.plot(xT_coordinates, yF1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xT_coordinates, yF2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(time_span, F_model, label = 'Fettreserv', linestyle="-", linewidth=lw, color=cb_palette1[5]) # Lägga till modellen
plt.legend((line3, line2, line1), ("Modell", "Högsta gräns", "Lägsta gräns"))
plt.xlabel("time", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Fettreserver")


# Sparar figur i plot constrains, glukos i muskeln
# Write the result to file
path_result_dir = "optimering/Bilder/test_plot_week17"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_fettreserve.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

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
import scipy.stats

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

    k1, k2, k3, k4, k5 = b
     
    # Says that the concentrations can't be lower than zero 
    x[x < 0] = 0.0 

    # Concentrations in the model as input 

    G, I, C, M, H = x 

    L= 5000 # Startvärde glukos i levern

    # Glucose plasma [1]
    dG = k4*C*I + k1*H - k2*G*I

    # Insulin plasma [2]
    dI = k3*G - k2*G*I

    # GLucose liver [3]
    dC = -k4*C*I + L

    # Glucose musle [4]
    dM = k2*G*I - k5*M

    # Glucose intake [5]
    dH = -k1*H

    return [dG, dI, dC, dM, dH]

def cost_function(b, yG_vec, yI_vec):
    if(any(b <= 0)):
        raise ValueError(f"{b}")

   # Calculates the target function for a model based on maximumlikelihood 

    # Start concentration, timespan   
    x0 = [30, 2.2e-8, 100, 60, 200]  # G, I, C, M, H, 
    
    #Injection of insulin
    inj = 7.3125E-07

    # Solve ODE-system until 20 minutes
    first_sol_G = integrate.solve_ivp(open_loop, time_span1, x0, method="Radau", args=(b, ), t_eval=tG_1)
    first_sol_I = integrate.solve_ivp(open_loop, time_span1, x0, method="Radau", args=(b, ), t_eval=tI_1) 
    
   # Simulate injection of insulin
    x1_G = [first_sol_G.y[0,-1], inj, first_sol_G.y[2,-1], first_sol_G.y[3,-1], first_sol_G.y[4,-1]]
    x1_I = [first_sol_I.y[0,-1], inj, first_sol_I.y[2,-1], first_sol_I.y[3,-1], first_sol_I.y[4,-1]]

    # Solve ODE-system after injection
    second_sol_G = integrate.solve_ivp(open_loop, time_span_G2, x1_G, method="Radau", args=(b, ), t_eval = tG_2)
    second_sol_I = integrate.solve_ivp(open_loop, time_span_I2, x1_I, method="Radau", args=(b, ), t_eval = tI_2)

    # The solution for the ODE-system over tG_vec and tI_vec
    sol_G = np.concatenate([first_sol_G.y, second_sol_G.y], axis = 1)
    sol_I = np.concatenate([first_sol_I.y, second_sol_I.y], axis = 1)
     
    # Solve ODE-system qualitative
    first_sol_qual = integrate.solve_ivp(open_loop, [0,20], x0, method="Radau", args=(b, ))

    # Simulate the injection
    x2 = [first_sol_qual.y[0, -1], inj, first_sol_qual.y[2, -1], first_sol_qual.y[3, -1], first_sol_qual.y[4, -1]]

    # Solve ODE-system after 20 miunutes with injection
    second_sol_qual = integrate.solve_ivp(open_loop, [20,240], x2, method = "Radau", args = (b, ))

    sol_qual = np.concatenate([first_sol_qual.y, second_sol_qual.y], axis = 1)

    G_model = sol_qual[0]
    I_model = sol_qual[1]
    C_model = sol_qual[2]
    M_model = sol_qual[3]
    H_model = sol_qual[4]

    # Extract G and I model concentrations at t-points tG_vec and tI_vec
    yG_model = sol_G[0] 
    yI_model = sol_I[1] 

    # Build bounds for the concentrations and punnish the cost-func. if it cross the bounds
    squared_sum = 0.0

    range_G = [0, 500] # mM 
    range_I = [0, 1.42e-6] #pM 
    range_C = [0, 10000] # mmol 
    range_M = [0, 500] # mmol
    range_H = [0, 500] # mmol

    penalty = 10000

    if any(G_model > np.max(range_G)):
       squared_sum += penalty
    if any(G_model < np.min(range_G)):
       squared_sum += penalty
    if any(I_model > np.max(range_I)):
       squared_sum += penalty
    if any(I_model < np.min(range_I)):
       squared_sum += penalty
    if any(C_model > np.max(range_C)):
       squared_sum += penalty
    if any(C_model < np.min(range_C)):
       squared_sum += penalty
    if any(M_model > np.max(range_M)):
       squared_sum += penalty
    if any(M_model < np.min(range_M)):
       squared_sum += penalty
    if any(H_model > np.max(range_H)):
       squared_sum += penalty
    if any(H_model < np.min(range_H)):
       squared_sum += penalty
    

    # Calculate cost-function  
    squared_sum = np.sum((yG_model - yG_vec)**2) +np.sum((yI_model -  yI_vec)**2)

    return squared_sum 
 
n_samples = 1000

n_a = -3
n_b = 3

a = 10 ** n_a
b = 10 ** n_b

y = scipy.stats.loguniform.rvs(a, b, size=(5, n_samples))

minimum = (np.inf, None)

for n in tqdm(range(n_samples)):

    cost = cost_function(y[:,n], cG_vec, cI_vec)

    if cost < minimum[0]:
        minimum =(cost, y[:,n])

print("Minimal costfunction and its parameters Log")
print(minimum)

# Start concentration, timespan   
x0 = [30, 2.2e-8, 100, 60, 200]  # G, I, C, M, H 
time_span_G = [tG_vec[0], tG_vec[-1]] 
time_span_I = [tI_vec[0], tI_vec[-1]]

#Injection
inj = 7.3125E-07

# Solve ODE-system qualitative
first_sol_qual = integrate.solve_ivp(open_loop, [0,20], x0, method="Radau", args=(minimum[1], ))

# Simulate the injection
x2 = [first_sol_qual.y[0, -1], inj,  first_sol_qual.y[2, -1], first_sol_qual.y[3, -1], first_sol_qual.y[4, -1]]

# Solve ODE-system after 20 miunutes with injection
second_sol_qual = integrate.solve_ivp(open_loop, [20,tG_vec[-1]], x2, method = "Radau", args = (minimum[1], ))

sol_qual = np.concatenate([first_sol_qual.y, second_sol_qual.y], axis = 1)
time_span = np.concatenate([first_sol_qual.t, second_sol_qual.t])

G_model = sol_qual[0]
I_model = sol_qual[1]
C_model = sol_qual[2]
M_model = sol_qual[3]
H_model = sol_qual[4]

# Time span
xT_coordinates = [tG_vec[0],tG_vec[-1]]

# Constrains glukos (G)
yG1_coordinates = [0,0] #mM (human)
yG2_coordinates = [50,50]   #mM (human)

# Constrains insulin (I)
yI1_coordinates = [0,0]   #pM (human)
yI2_coordinates = [1.42e-6,1.42e-6] #pM (human)


# plotta glukos
lw = 2.0
plot1 = plt.figure(1)
# G_plot = plt.subplot(121)
line1, = plt.plot(xT_coordinates, yG1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xT_coordinates, yG2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(data_G['time'].values, data_G['conc'].values, label = 'Glukos', linestyle="-", linewidth=lw, color=cb_palette1[7])
line4, = plt.plot(time_span, G_model, label = 'Glukos', linestyle="-", linewidth=lw, color=cb_palette1[5])
plt.legend((line4, line3, line2, line1), ("Modell", "Data", "Högsta gräns","Lägsta gräns"))
plt.xlabel("Tid [min]", fontsize=12), plt.ylabel("Konc.[mM]", fontsize=12)
plt.title("Glukos i plasma")

# # Residual plot for glucose
# G_res = plt.subplot(122)
# difference_G = cG_vec - G_model
# plt.scatter(difference_G, G_model, s = 10 , color = cb_palette1[1])

# Sparar figur i plot constrains, glukos
# Write the result to file
path_result_dir = "optimering/Bilder/no_op_plot_week16"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.makedirs(path_result_dir, exist_ok=True)  # Create a new directory if not existing
path_fig = path_result_dir + "/no_op_plot_glucose.pdf"
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
plt.xlabel("Tid [min]", fontsize=12), plt.ylabel("Konc. [mM]", fontsize=12)
plt.title("Insulin i plasma")

# # Residual plot for insulin
# I_res = plt.subplot(122)
# difference_I = cI_vec - I_model
# plt.scatter(difference_I, I_model, s = 10 , color = cb_palette1[1])

# Sparar figur i plot constrains, insulin
# Write the result to file
path_result_dir = "optimering/Bilder/no_op_plot_week16"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/no_op_plot_insulin.pdf"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta Glukos i lever 
lw = 2.0
plot1 = plt.figure(3)
# line1, = plt.plot(xT_coordinates, yC1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
# line2, = plt.plot(xT_coordinates, yC2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(time_span, C_model, label = 'Modell', linestyle="-", linewidth=lw, color=cb_palette1[5]) # Lägga till modellen
plt.legend()
plt.xlabel("Tid [min]", fontsize=12), plt.ylabel("Konc. [mM]", fontsize=12)
plt.title("Glukos i levern")

# Sparar figur i plot constrains, glukos i levern
# Write the result to file
path_result_dir = "optimering/Bilder/no_op_plot_week16"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/no_op_plot_glukoslevern.pdf"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta Glukos i muskeln 
lw = 2.0
plot1 = plt.figure(4)
# line1, = plt.plot(xT_coordinates, yM1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
# line2, = plt.plot(xT_coordinates, yM2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(time_span, M_model, label = 'Modell', linestyle="-", linewidth=lw, color=cb_palette1[5]) # Lägga till modellen
plt.legend()
plt.xlabel("Tid [min]", fontsize=12), plt.ylabel("Konc. [mM]", fontsize=12)
plt.title("Glukos i muskeln")


# Sparar figur i plot constrains, glukos i muskeln
# Write the result to file
path_result_dir = "optimering/Bilder/no_op_plot_week16"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/no_op_plot_glukosmuskeln.pdf"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)


# plotta Glukos inktake 
lw = 2.0
plot1 = plt.figure(5)
# line1, = plt.plot(xT_coordinates, yH1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
# line2, = plt.plot(xT_coordinates, yH2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(time_span, H_model, label = 'Modell', linestyle="-", linewidth=lw, color=cb_palette1[5]) # Lägga till modellen
plt.legend()
plt.xlabel("Tid [min]", fontsize=12), plt.ylabel("Konc. [mM]", fontsize=12)
plt.title("Glukos intag")


# Sparar figur i plot constrains, glukos i muskeln
# Write the result to file
path_result_dir = "optimering/Bilder/no_op_plot_week16"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/no_op_plot_glukos_intag.pdf"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)
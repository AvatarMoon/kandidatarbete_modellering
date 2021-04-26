import numpy as np 
import scipy.integrate as integrate 
from scipy.optimize import minimize
from scipy.optimize import Bounds
import matplotlib.pyplot as plt 
import math 
import pandas as pd
import os
import lhsmdu

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


def open_loop(t,x,b): 
    k1, k2, k3, k4, k5 = b
     
    # Says that the concentrations can't be lower than zero 
    x[x < 0] = 0.0 

    # Concentrations in the model as input 
    G, I, C, M = x 

     # Startvärden
    H= 600 # Intag glukos 
    L= 500 # Startvärde glukos i levern

    # Glucose plasma [1]
    dG = k4*I*C + k1*H - k2*M

    # Insulin plasma [2]
    dI = -k4*I*C + k3*G

    # GLucose liver [4]
    dC = -k4*I*C + L

    # Glucose musle [5]
    dM = k2*G - k5*M

    return [dG, dI, dC, dM]

 

def cost_function(b, yG_vec, yI_vec): 

    # Start concentration, timespan   
    x0 = [60, 5, 34, 3]  # G, I, C, M 
    time_span_G = [tG_vec[0], tG_vec[-1]] 
    time_span_I = [tI_vec[0], tI_vec[-1]] 
    

    # Step 1: Solve ODE-system at points tG_vec
    sol_G = integrate.solve_ivp(open_loop, time_span_G, x0, method="LSODA", args=(b, ), t_eval=tG_vec) 
    sol_I = integrate.solve_ivp(open_loop, time_span_I, x0, method="LSODA", args=(b, ), t_eval=tI_vec) 

    
    # step 2: Solve ODE-system qualitative
    sol_qual = integrate.solve_ivp(open_loop, time_span_G, x0, method="LSODA", args=(b, ))

    G_model = sol_qual.y[0]
    I_model = sol_qual.y[1]
    C_model = sol_qual.y[2]
    M_model = sol_qual.y[3]

    # Step 3: Extract G and I model concentrations at t-points tG_vec and tI_vec
    yG_model = sol_G.y[0] 
    yI_model = sol_I.y[1] 

    # Step 4 : Build bounds for the concentrations and punnish the cost-func. if they go cross the bounds
    squared_sum = 0.0

    range_G = [0, 50] # mM 
    range_I = [0, 5000] #pM 
    range_C = [0, 100] # mmol 
    range_M = [0, 140] # mmol 

    if any(G_model) > np.max(range_G):
       squared_sum += 100
    if any(G_model) < np.min(range_G):
        squared_sum += 100
    if any(I_model) > np.max(range_I):
       squared_sum += 100
    if any(I_model) < np.min(range_I):
        squared_sum += 100
    if any(C_model) > np.max(range_C):
       squared_sum += 100
    if any(C_model) < np.min(range_C):
        squared_sum += 100
    if any(M_model) > np.max(range_M):
       squared_sum += 100
    if any(M_model) < np.min(range_M):
        squared_sum += 100
    

    # Step 5: Calculate cost-function  
    squared_sum = np.sum((yG_model - yG_vec))**2+np.sum((yI_model -  yI_vec)**2) 

    return squared_sum 

## Hypercube set up
randSeed = 2 # random number of choice
lhsmdu.setRandomSeed(randSeed) # Latin Hypercube Sampling with multi-dimensional uniformity
start = np.array(lhsmdu.sample(5, 4)) # Latin Hypercube Sampling with multi-dimensional uniformity (parameters, samples)

para, samples = start.shape

## intervals for the parameters
para_int = [0, 500]

minimum = (np.inf, None)

# Bounds for the model
bound_low = np.array([0, 0, 0, 0, 0])
bound_upp = np.repeat(np.inf, para)
bounds = Bounds(bound_low, bound_upp)


for n in range(samples):
    k1 = start[0,n] * para_int[1]
    k2 = start[1,n] * para_int[1]
    k3 = start[2,n] * para_int[1]
    k4 = start[3,n] * para_int[1]
    k5 = start[4,n] * para_int[1]
    
    res = minimize(cost_function, [k1, k2, k3, k4, k5], method='Powell', args = (cG_vec, cI_vec), bounds=bounds) #lägg in constraints här 

    if res.fun < minimum[0]:
        minimum = (res.fun, res.x)

# Hämta modellen
    # Start concentration, timespan   
    x0 = [60, 5, 34, 3]  # G, I, C, M 
    time_span_G = [tG_vec[0], tG_vec[-1]] 
    time_span_I = [tI_vec[0], tI_vec[-1]] 
    

    # Step 1: Solve ODE-system at points tG_vec
    sol_G = integrate.solve_ivp(open_loop, time_span_G, x0, method="LSODA", args=(minimum[1], ), t_eval=tG_vec) 
    sol_I = integrate.solve_ivp(open_loop, time_span_I, x0, method="LSODA", args=(minimum[1], ), t_eval=tI_vec) 

    sol_qual = integrate.solve_ivp(open_loop, time_span_G, x0, method="LSODA", args=(minimum[1], ))

    G_model = sol_qual.y[0]
    I_model = sol_qual.y[1]
    C_model = sol_qual.y[2]
    M_model = sol_qual.y[3]

# Print some statistics  
print("Optimal value found via Powells-method:") 
print(minimum[1]) 
print("Value of cost-function") 
print(minimum[0]) 

# Plot model, data and constrains

# Constrains glukos (G)
xG_coordinates = [tG_vec[0],tG_vec[-1]]
yG1_coordinates = [0,0] #mM (human)
yG2_coordinates = [50,50]   #mM (human)

# Constrains insulin (I)
xI_coordinates = [tI_vec[0],tI_vec[-1]]
yI1_coordinates = [0,0]   #pM (human)
yI2_coordinates = [5000,5000] #pM (human)
 
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
path_result_dir = "optimering/Bilder/plot_week16_open_loop_model"
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
path_result_dir = "optimering/Bilder/plot_week16_open_loop_model"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_insulin.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta Glukos i lever 
lw = 2.0
plot1 = plt.figure(3)
line1, = plt.plot(xC_coordinates, yC1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xC_coordinates, yC2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(data_G['time'].values, C_model, label = 'Glukos i levern', linestyle="-", linewidth=lw, color=cb_palette1[5]) # Lägga till modellen
plt.legend((line3, line2, line1), ("Modell", "Högsta gräns", "Lägsta gräns"))
plt.xlabel("time", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Glukos i levern")

# Sparar figur i plot constrains, glukos i levern
# Write the result to file
path_result_dir = "optimering/Bilder/plot_week16_open_loop_model"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_glukoslevern.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta Glukos i muskeln 
lw = 2.0
plot1 = plt.figure(4)
line1, = plt.plot(xM_coordinates, yM1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xM_coordinates, yM2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(data_G['time'].values, M_model, label = 'Glukos i muskeln', linestyle="-", linewidth=lw, color=cb_palette1[5]) # Lägga till modellen
plt.legend((line3, line2, line1), ("Modell", "Högsta gräns", "Lägsta gräns"))
plt.xlabel("time", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Glukos i muskeln")
plt.show()

# Sparar figur i plot constrains, glukos i muskeln
# Write the result to file
path_result_dir = "optimering/Bilder/plot_week16_open_loop_model"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_glukosmuskeln.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

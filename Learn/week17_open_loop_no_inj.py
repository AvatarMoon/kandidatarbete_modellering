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

def open_loop(t,x,b): 

    k1, k2, k3, k4, k5, k6, k7, k8 = b
     
    # Says that the concentrations can't be lower than zero 
    x[x < 0] = 0.0 

    # Concentrations in the model as input 

    G, I, C, M, H, E, F = x 

    #L= 5000 # Startvärde glukos i levern


    # Glucose plasma [1]
    dG = k4*C + k1*H - k2*G*I

    # Insulin plasma [2]
    dI = k3*G - k2*G*I

    # GLucose liver [3]
    dC = -k4*C + k6*E + k7*F 

    # Glucose musle [4]
    dM = k2*G*I - k5*M

    # Glucose intake [5]
    dH = -k1*H

    # Glucagon [6]
    dE = k8 - k2*G*E

    # Fettreserv [7]
    dF = -k7*F

    return [dG, dI, dC, dM, dH, dE, dF]


def cost_function(b, yG_vec, yI_vec):
    if(any(b <= 0)):
        raise ValueError(f"{b}")

   # Calculates the target function for a model based on maximumlikelihood 

    # Start concentration, timespan   
    x0 = [30, 2.2e-8, 34, 60, 70, 50, 400]  # G, I, C, M, H, E, F 
    time_span_G = [tG_vec[0], tG_vec[-1]] 
    time_span_I = [tI_vec[0], tI_vec[-1]] 


    # Step 1: Solve ODE-system at points tG_vec
    sol_G = integrate.solve_ivp(open_loop, time_span_G, x0, method="Radau", args=(b, ), t_eval=tG_vec) 
    sol_I = integrate.solve_ivp(open_loop, time_span_I, x0, method="Radau", args=(b, ), t_eval=tI_vec) 

    
    # step 2: Solve ODE-system qualitative
    sol_qual = integrate.solve_ivp(open_loop, time_span_G, x0, method="Radau", args=(b, ))

    G_model = sol_qual.y[0]
    I_model = sol_qual.y[1]
    C_model = sol_qual.y[2]
    M_model = sol_qual.y[3]
    H_model = sol_qual.y[4]
    E_model = sol_qual.y[5]
    F_model = sol_qual.y[6]


    # Extract G and I model concentrations at t-points tG_vec and tI_vec
    yG_model = sol_G.y[0] 
    yI_model = sol_I.y[1] 

    # Build bounds for the concentrations and punnish the cost-func. if it cross the bounds
    squared_sum = 0.0

    range_G = [0, 500] # mM 
    range_I = [0, 1.4e-6] #pM 
    range_C = [0, 1000] # mmol 
    range_M = [0, 1000] # mmol
    range_H = [0, 500] # mmol
    range_E = [0, 500] # mmol
    range_F = [0, 500] # mmol

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
    if any(E_model > np.max(range_E)):
       squared_sum += penalty
    if any(E_model < np.min(range_E)):
       squared_sum += penalty
    if any(F_model > np.max(range_F)):
       squared_sum += penalty
    if any(F_model < np.min(range_F)):
       squared_sum += penalty
    

    # Calculate cost-function  
    squared_sum += np.sum((yG_model - yG_vec)**2) + np.sum((yI_model -  yI_vec)**2) 

    return squared_sum 


## Hypercube set up
randSeed = 2 # random number of choice
lhsmdu.setRandomSeed(randSeed) # Latin Hypercube Sampling with multi-dimensional uniformity
start = np.array(lhsmdu.sample(8, 10)) # Latin Hypercube Sampling with multi-dimensional uniformity (parameters, samples)

para, samples = start.shape

## intervals for the parameters
para_int = [0, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000]
#           0     1      2      3    4   5   6   7     8
minimum = (np.inf, None)

# Bounds for the model

bound_low = np.array([0, 0, 0, 0, 0, 0, 0, 0])

bound_upp = np.repeat(np.inf, para)
bounds = Bounds(bound_low, bound_upp)

fig = plt.figure() 

os.makedirs("logs17", exist_ok=True)
filename = f"logs17/{datetime.datetime.utcnow()}.log"
filename = filename.replace(" ","_")
filename = filename.replace(":",".")

for n in tqdm(range(samples)):

    k1 = start[0,n] * para_int[2]
    k2 = start[1,n] * para_int[2]
    k3 = start[2,n] * para_int[2]
    k4 = start[3,n] * para_int[2]
    k5 = start[4,n] * para_int[2]
    k6 = start[5,n] * para_int[2]
    k7 = start[6,n] * para_int[2]
    k8 = start[7,n] * para_int[2]
    
    try:      
        res = minimize(
            fun = cost_function, 
            x0 = [k1, k2, k3, k4, k5, k6, k7, k8], 
            method ='Powell', 
            args = (cG_vec, cI_vec), 
            bounds = bounds,
            tol = 0.1,
            options = {'disp' : True},
            ) #lägg in constraints här 
    except ValueError as err:
        msg = f"Start values: {[k1, k2, k3, k4, k5, k6, k7, k8]}\nLed to negative solution: {err}"
        with open(filename, "a") as f:
            f.writelines([f"Failed!\n", msg])
        print(msg)
        continue

    if res.fun < minimum[0]:
        minimum = (res.fun, res.x)
    with open(filename, "a") as f:
        f.writelines([f"Success!", f"Found solution: {res.x}"])


# Hämta modellen
# Start concentration, timespan   

x0 = [30, 2.2e-8, 34, 60, 70, 50, 400]  # G, I, C, M, H, E, F 

time_span_G = [tG_vec[0], tG_vec[-1]] 
time_span_I = [tI_vec[0], tI_vec[-1]] 

# Solve ODE-system qualitative
sol_qual = integrate.solve_ivp(open_loop, time_span_G, x0, method="Radau", args=(minimum[1], ))


G_model = sol_qual.y[0]
I_model = sol_qual.y[1]
C_model = sol_qual.y[2]
M_model = sol_qual.y[3]
H_model = sol_qual.y[4]
E_model = sol_qual.y[5]
F_model = sol_qual.y[6]


# Print some statistics  
print("Optimal value found via Powells-method:") 
print(minimum[1]) 
print("Value of cost-function") 
print(minimum[0]) 

def sensitivity(b,t,x):
   # Calcualte the sensitivity matrix using the optimal 
    h = np.sqrt(np.finfo(float).eps) # Maskintoleransen, vår steglängd för finita differen 
    b_par = len(b)
    t_len = len(t)
    # Sensitivity analysis for each time step
    S = np.zeros([b_par, t_len * len(x)])
    time_span = [t[0], t[-1]]

    for n in range(len(b)):
        b1 = b.copy() 
        b2 = b.copy()  
        b1[n] += h 
        b2[n] -= h

        Sol_high = integrate.solve_ivp(open_loop, time_span, x, method='LSODA', args=(b1, ), t_eval = t)
        Sol_low = integrate.solve_ivp(open_loop, time_span, x, method='LSODA', args=(b2, ), t_eval= t)
        
        Sol_diff = (Sol_high.y-Sol_low.y)/(2*h)

        S[n,:] = Sol_diff.reshape(t_len*len(x))

    return S

S = sensitivity(minimum[1], tG_vec , x0)

# Fisher matrix to make the covariance matrix
Fisher = 2 * S @ S.transpose()

cov_mat = Fisher.transpose()

# Identification of the parameters
d_cov = np.diag(cov_mat)

var_coeff = np.square(d_cov)/minimum[1]

print('Identification for each parameters')
print(var_coeff)

### ~~~~~~ Plot model, data, constrains and residual ~~~~~~~ ###

# Time span
time_span = np.linspace(tG_vec[0], tG_vec[-1], len(G_model))
xT_coordinates = [tG_vec[0],tG_vec[-1]]

# Constrains glukos (G)
yG1_coordinates = [0,0] #mM (human)
yG2_coordinates = [50,50]   #mM (human)

# Constrains insulin (I)
yI1_coordinates = [0,0]   #pM (human)
yI2_coordinates = [4000,4000] #pM (human)
 
# Constrains glukos i lever (C) 
yC1_coordinates = [0,0]  # mmol (human)
yC2_coordinates = [1000,1000]  # mmol (human)

# Constrains glukos i muskel (M) 
yM1_coordinates = [0,0]    # mmol (human)
yM2_coordinates = [1400,1400]  # mmol (human)

# Constrains glucose intake (H)
yH1_coordinates = [0,0]    
yH2_coordinates = [200,200]  

# Constrains glucagon plasma (E)
yE1_coordinates = [0,0]    
yE2_coordinates = [75, 75]  

# Constrains Fettreserver (F)
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
plt.xlabel("tid (min)", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Glukos i plasma")


# Write the result to file
path_result_dir = "optimering/Bilder/plot_week17_no_inj_model"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
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
plt.xlabel("tid", fontsize=12), plt.ylabel("Insulin koncentration", fontsize=12)
plt.title("Insulin i plasma")


# Sparar figur i plot constrains, insulin
# Write the result to file
path_result_dir = "optimering/Bilder/plot_week17_no_inj_model"
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
plt.xlabel("tid", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Glukos i levern")

# Sparar figur i plot constrains, glukos i levern
# Write the result to file
path_result_dir = "optimering/Bilder/plot_week17_no_inj_model"
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
plt.xlabel("tid", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Glukos i muskeln")


# Sparar figur i plot constrains, glukos i muskeln
# Write the result to file
path_result_dir = "optimering/Bilder/plot_week17_no_inj_model"
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
plt.xlabel("tid", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Glukos intag")


# Sparar figur i plot constrains, glukos i muskeln
# Write the result to file
path_result_dir = "optimering/Bilder/plot_week17_no_inj_model"
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
path_result_dir = "optimering/Bilder/plot_week17_no_inj_model"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_glukagon_i_plasma.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta Fettreserver
lw = 2.0
plot1 = plt.figure(7)
line1, = plt.plot(xT_coordinates, yF1_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[1])
line2, = plt.plot(xT_coordinates, yF2_coordinates, linestyle=":", linewidth=lw, color=cb_palette1[3])
line3, = plt.plot(time_span, F_model, label = 'Fettreserver', linestyle="-", linewidth=lw, color=cb_palette1[5]) # Lägga till modellen
plt.legend((line3, line2, line1), ("Modell", "Högsta gräns", "Lägsta gräns"))
plt.xlabel("tid", fontsize=12), plt.ylabel("Glukos koncentration", fontsize=12)
plt.title("Fettreserver")


# Sparar figur i plot constrains, glukos i muskeln
# Write the result to file
path_result_dir = "optimering/Bilder/plot_week17_no_inj_model"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_Fettreserver.jpg"
path_fig = path_result_dir + "/plot_glucagon_plasma.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)
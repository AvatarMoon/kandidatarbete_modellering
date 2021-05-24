import lhsmdu
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd


cb_palette1 = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

# # ~~~ Plot data

# # get data 
# data_G = pd.read_csv ("data_horses/glukos_FF_training.csv", sep=';')
# data_I = pd.read_csv ("data_horses/insulin_FF_training.csv", sep=';')
# data_G = data_G.sort_values(by=['time'])
# data_I = data_I.sort_values(by=['time'])


# # Extract times- and concentration-vectors
# tG_vec = data_G['time'].values
# tI_vec = data_I['time'].values  
# cG_vec = data_G['conc'].values
# cI_vec = data_I['conc'].values 

# # plotta Glukos 
# lw = 2.0
# plot1 = plt.figure(1)
# line1, = plt.plot(tG_vec, cG_vec, linestyle="-", linewidth=lw, color=cb_palette1[7]) # Lägga till modellen
# plt.xlabel("Tid [min]", fontsize=15), plt.ylabel("Konc. [mM]", fontsize=15)
# plt.title("Data glukos i plasma", fontsize=15)

# # Sparar figur i plot constrains, glukos i muskeln
# # Write the result to file
# path_result_dir = "data_horses/bilder"
# # Check if directory exists
# if not os.path.isdir(path_result_dir):
#     os.mkdir(path_result_dir)  # Create a new directory if not existing
# path_fig = path_result_dir + "/plot_glukos.pdf"
# print("path_fig = {}".format(path_fig))
# plt.savefig(path_fig)

# # plotta Inslulin
# lw = 2.0
# plot1 = plt.figure(2)
# line1, = plt.plot(tI_vec, cI_vec, linestyle="-", linewidth=lw, color=cb_palette1[7]) # Lägga till modellen
# plt.xlabel("Tid [min]", fontsize=15), plt.ylabel("Konc. [mM]", fontsize=15)
# plt.title("Data insulin i plasma", fontsize=15)

# # Sparar figur i plot constrains, glukos i muskeln
# # Write the result to file
# path_result_dir = "data_horses/bilder"
# # Check if directory exists
# if not os.path.isdir(path_result_dir):
#     os.mkdir(path_result_dir)  # Create a new directory if not existing
# path_fig = path_result_dir + "/plot_insulin.pdf"
# print("path_fig = {}".format(path_fig))
# plt.savefig(path_fig)


## ~~~ Plot random vs hypercube
randSeed = 2 # random number of choice
lhsmdu.setRandomSeed(randSeed) # Latin Hypercube Sampling with multi-dimensional uniformity
l = lhsmdu.sample(2,50) # Latin Hypercube Sampling of two variables, and 10 samples each.
k = lhsmdu.createRandomStandardUniformMatrix(2,50) # Monte Carlo Sampling

for i in range(1,50//5+1):
    ## Hypercube set up

    fig1 = plt.figure()
    plt.axis([0, 1, 0, 1])
    plt.scatter([k[0,0:i*5]], [k[1, 0:i*5]], color=cb_palette1[0])
    plt.title(f" {i*5} Likformigt slumpad samplingar", fontsize = '15')

    # Write the result to file
    path_result_dir = "Modeller/Bilder"
    # Check if directory exists
    if not os.path.isdir(path_result_dir):
        os.mkdir(path_result_dir)  # Create a new directory if not existing
    path_fig = path_result_dir + f"/random_{i*5}.pdf"
    print("path_fig = {}".format(path_fig))
    plt.savefig(path_fig)

    fig2 = plt.figure() 
    plt.axis([0, 1, 0, 1])
    plt.scatter([l[0,0:i*5]], [l[1,0:i*5]], color=cb_palette1[5])
    plt.title(f"{i*5} samplingar från LHS", fontsize = '15')

    # Write the result to file
    path_result_dir = "Modeller/Bilder"
    # Check if directory exists
    if not os.path.isdir(path_result_dir):
        os.mkdir(path_result_dir)  # Create a new directory if not existing
    path_fig = path_result_dir + f"/LHS_{i*5}.pdf"
    print("path_fig = {}".format(path_fig))
    plt.savefig(path_fig)

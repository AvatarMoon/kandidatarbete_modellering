import matplotlib.pyplot as plt
import numpy 
import pandas as pd
import os

# Colour-blind friendly palette (use nice colors)
cb_palette1 = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

# [Mörkblå, Gul, Ockraröd, Mörkgrön ,Olivgrön ,Ljusbeige (dålig), Ljusgult (dålig), Gul ]
cb_palette2 = ["#F4E3AF", "#F1CB6F", "#E16F3B", "#2D4E63", "#899964", "#F4E3AF", "#F1CB6F", "#E16F3B"]

# Hämta data
data_G = pd.read_csv ('data_horses\Glukos_new_FFaraber.csv', sep=';')
data_I = pd.read_csv ('data_horses\insulin_new_FFaraber.csv', sep=';')
data_G = data_G.sort_values(by = ['tid'])
data_I = data_I.sort_values(by = ['tid'])

# Tidsvektorer
tG_vec = data_G['tid'].values
tI_vec = data_I['tid'].values

# Constrains glukos
xG_coordinates = [tG_vec[0],tG_vec[-1]]
yG1_coordinates = [4.5,4.5] #mM (human)
yG2_coordinates = [11,11]   #mM (human)

# Constrains insulin
xI_coordinates = [tI_vec[0],tI_vec[-1]]
yI1_coordinates = [38,38]   #pM (human)
yI2_coordinates = [400,400] #pM (human)

# plotta glukos
lw = 2.0
plot1 = plt.figure(1)
line1, = plt.plot(xG_coordinates, yG1_coordinates, linestyle=":", linewidth=lw, color=cb_palette2[1])

plot1 = plt.figure(1)
line2, = plt.plot(xG_coordinates, yG2_coordinates, linestyle=":", linewidth=lw, color=cb_palette2[3])

plot1 = plt.figure(1)
line3, = plt.plot(data_G['tid'].values, data_G['konc'].values, label = 'Glukos', linestyle="-", linewidth=lw, color=cb_palette2[2])
plt.legend((line1, line2, line3), ("Lägsta gräns", "Högsta gräns", "Data"))
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
line1, = plt.plot(xI_coordinates, yI1_coordinates, linestyle=":", linewidth=lw, color=cb_palette2[1])

plot1 = plt.figure(2)
line2, = plt.plot(xI_coordinates, yI2_coordinates, linestyle=":", linewidth=lw, color=cb_palette2[3])

plot1 = plt.figure(2)
line3, = plt.plot(data_I['tid'].values, data_I['conc'].values, label = 'Insulin', linestyle="-", linewidth=lw, color=cb_palette2[4])
plt.legend((line1, line2, line3), ("Lägsta gräns", "Högsta gräns", "Data"))
plt.xlabel("time", fontsize=12), plt.ylabel("Insulin koncentration", fontsize=12)
plt.title("Insulin in plasma")
plt.show()

# Sparar figur i plot constrains, insulin
# Write the result to file
path_result_dir = "optimering/Bilder/plot_constrains"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_insulin.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)
import matplotlib.pyplot as plt
import numpy 
import pandas as pd
import os

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
plot1 = plt.figure(1)
plt.plot(xG_coordinates, yG1_coordinates)
plot1 = plt.figure(1)
plt.plot(xG_coordinates, yG2_coordinates)
plot1 = plt.figure(1)
plt.plot(data_G['tid'].values, data_G['konc'].values, label = 'Glukos')
plt.legend()

# Sparar figur i plot constrains
# Write the result to file
path_result_dir = "optimering/Bilder/plot_constrains"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_glucose.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

# plotta insulin
plot1 = plt.figure(2)
plt.plot(xI_coordinates, yI1_coordinates)
plot1 = plt.figure(2)
plt.plot(xI_coordinates, yI2_coordinates)
plot1 = plt.figure(2)
plt.plot(data_I['tid'].values, data_I['conc'].values, label = 'Insulin')
plt.legend()

# Sparar figur i plot constrains
# Write the result to file
path_result_dir = "optimering/Bilder/plot_constrains"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/plot_insulin.jpg"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)
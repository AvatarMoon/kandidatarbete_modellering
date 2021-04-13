import matplotlib.pyplot as plt
import pandas as pd

data_G = pd.read_csv ('data_horses\Glukos_new_FFaraber.csv', sep=';')
data_I = pd.read_csv ('data_horses\insulin_new_FFaraber.csv', sep=';')

tG_vec = data_G['tid'].values
tI_vec = data_I['tid'].values

# Constrains glukos
xG_coordinates = [tG_vec[0],tG_vec[-1]]
yG1_coordinates = [4.5,4.5] #mM
yG2_coordinates = [11,11]   #mM

# Constrains insulin
xI_coordinates = [tI_vec[0],tI_vec[-1]]
yI1_coordinates = [38,38]   #pM
yI2_coordinates = [400,400] #pM

# plotta glukos
plot1 = plt.figure(1)
plt.plot(xG_coordinates, yG1_coordinates)
plot1 = plt.figure(1)
plt.plot(xG_coordinates, yG2_coordinates)
plot1 = plt.figure(1)
plt.plot(data_G['tid'].values, data_G['konc'].values, label = 'Glukos')
plt.legend()
plt.show()

# plotta insulin
plot1 = plt.figure(2)
plt.plot(xI_coordinates, yI1_coordinates)
plot1 = plt.figure(2)
plt.plot(xI_coordinates, yI2_coordinates)
plot1 = plt.figure(2)
plt.plot(data_I['tid'].values, data_I['conc'].values, label = 'Insulin')
plt.legend()
plt.show()
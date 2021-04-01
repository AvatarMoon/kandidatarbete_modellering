import matplotlib.pyplot as plt
import pandas as pd

data_G = pd.read_csv ("data_horses/Glukos_FFaraber.csv", sep=';')
data_I = pd.read_csv ("data_horses/Insulin_FFaraber.csv", sep=';')

plt.plot(data_G['tid'].values, data_G['konc'].values, label = 'Glukos')
plt.plot(data_I['tid'].values, data_I['konc'].values, label = 'Insulin')
plt.legend()
plt.show()
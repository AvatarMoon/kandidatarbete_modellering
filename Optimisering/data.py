import numpy as np
import os
import csv

data_path = os.listdir('data_horses')

def read_csv(data_path):
    data = []
    with open(data_path) as file:
        for line in file.readlines()[1:]:
            entries = line.strip().split(";")
            value = entries[1]
            data.append(value)
    return data

conc = []
for x in data_path:
    values_list = []
    data = read_csv('data_horses/' + x)
    conc.append(data)
    

g_data = conc[0]
i_data = conc[1]

print(g_data)
print(i_data)

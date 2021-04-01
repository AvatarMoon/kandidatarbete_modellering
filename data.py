import numpy as np
import os
import csv

data_path = os.listdir('data_horses')

conc = []
for x in data_path:
    values_list = []
    with open('data_horses/' + x) as file:
        for line in file.readlines()[1:]:
            entries = line.strip().split(";")
            value = entries[1]
            values_list.append(value)
    conc.append(values_list)

g_data = conc[0]
i_data = conc[1]

print(g_data)
print(i_data)

import numpy as np
import os
import csv
import pandas as pd

data_path = os.listdir('data_horses')

conc = [] 

for x in data_path:
    # data = pd.read_csv('data_horses/' + x, sep = ';')
    # conc = data.to_string()
    values_list = []
    with open('data_horses/' + x) as file:
        for line in file.readlines()[1:]:
            entries = line.strip().split(";")
            value = entries[1]
            values_list.append(value)
    conc.append(values_list)

print(conc)

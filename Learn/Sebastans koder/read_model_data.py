import numpy as np

file_path = "mmc1/SimulationFiles/Detailed models/M1 hierarchical model/plotAllGoodValues/allGoodValues.dat"

# Manually loading woth Python
with open(file_path) as f:
    values_list = []
    for line in f.readlines():
        entries = line.split(" ")
        values = [float(val) for val in entries[:-1]]
        values_list.append(values)

arr_python = np.array(values_list)

# Loading with numpy
arr = np.loadtxt(file_path)

print(np.array_equal(arr_python, arr))

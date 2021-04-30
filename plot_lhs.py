import lhsmdu
import numpy as np
import matplotlib.pyplot as plt

## ~~~ Plot random vs hypercube

l = lhsmdu.sample(2,20) # Latin Hypercube Sampling of two variables, and 10 samples each.
k = lhsmdu.createRandomStandardUniformMatrix(2,20) # Monte Carlo Sampling


fig1 = plt.figure(1)
plt.scatter([k[0]], [k[1]], color=cb_palette1[0])
plt.title("Random")

# Write the result to file
path_result_dir = "Modeller/"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/random.pdf"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)

fig2 = plt.figure(2)
plt.scatter([l[0]], [l[1]], color=cb_palette1[1])
plt.title("LHS")

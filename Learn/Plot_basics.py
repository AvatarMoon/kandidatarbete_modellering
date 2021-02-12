import numpy as np
import matplotlib.pyplot as plt
import os   # For creating folders and interacting with the operating system

# Colour-blind friendly palette (use nice colors)
cb_palette = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

# Generate some numpy format data
x_data = np.linspace(0, 1, num=100)
y1 = np.exp(x_data)
y2 = np.exp(x_data ** 2)
y3 = np.exp(x_data * 2)

# Plot the data
lw = 2.0
line1, = plt.plot(x_data, y1, linestyle="-", linewidth=lw, color=cb_palette[1])
line2, = plt.plot(x_data, y2, linestyle=":", linewidth=lw, color=cb_palette[2])
line3, = plt.plot(x_data, y3, linestyle="--", linewidth=lw, color=cb_palette[3])
plt.legend((line1, line2, line3), ("$e^x$", "$e^{x^2}$", "$e^{2x}$"))
plt.xlabel("x", fontsize=12), plt.ylabel("y", fontsize=12)
plt.title("Example graph")

# Write the result to file
path_result_dir = "Users/mauritzkoohler/Downloads.pdf"
# Check if directory exists
if not os.path.isdir(path_result_dir):
    os.mkdir(path_result_dir)  # Create a new directory if not existing
path_fig = path_result_dir + "/Example_graph.pdf"
print("path_fig = {}".format(path_fig))
plt.savefig(path_fig)




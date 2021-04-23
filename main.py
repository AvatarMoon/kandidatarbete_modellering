from optimering.data_loaders import load_glukos, load_insulin
from optimering.models import open_loop

# get data from horses for open_loop
tG, cG = load_glukos()
tI, cI = load_insulin()

#  1) Optimize an open_loop model on its parameters. 
optimal_parameters = optimize_open_loop(open_loop)
#  2) Makes a diagnositic on the model against the data
#  3) Makes a sesitivity analysis and identification on the model 
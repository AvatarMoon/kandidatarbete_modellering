# Göra om data till en matris
import numpy as np 
import scipy.io
mat_274 = scipy.io.loadmat('expdata_GLU_time', mdict=None, appendmat=True)

data = np.array(mat_274)
print(data)
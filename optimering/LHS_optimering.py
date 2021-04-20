import lhsmdu
import numpy as np
import numpy
import matplotlib.pyplot as plt

k = lhsmdu.sample(3, 20) # Latin Hypercube Sampling with multi-dimensional uniformity
k = np.array(k)


fig = plt.figure()
ax = fig.gca()
ax.set_xticks(numpy.arange(0,1,0.1))
ax.set_yticks(numpy.arange(0,1,0.1))

plt.scatter(k[0], k[1], color='g', s = 30)
plt.grid()
plt.show()
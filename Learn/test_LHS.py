import lhsmdu
import matplotlib.pyplot as plt
import numpy

l = lhsmdu.sample(2,10) # Latin Hypercube Sampling of two variables, and 10 samples each.
k = lhsmdu.createRandomStandardUniformMatrix(2,10) # Monte Carlo Sampling

fig = plt.figure()
ax = fig.gca()
ax.set_xticks(numpy.arange(0,1,0.1))
ax.set_yticks(numpy.arange(0,1,0.1))
plt.scatter(k[0], k[1], color="b", label="LHS-MDU")
plt.scatter(l[0], l[1], color="r", label="MC")
plt.grid()
plt.show()
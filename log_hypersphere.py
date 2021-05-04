import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

# Uniform
def hypersphere(n_dimensions, n_samples):
    rho = np.random.uniform(0, 1, n_samples) ** (1 / n_dimensions)
    r = np.abs(np.random.normal(size=(n_dimensions, n_samples)))
    c = np.sqrt(np.sum(r ** 2, axis=0))
    return rho * r / c


# Logarithmic radius
def log_hypersphere(a, b, n_dimensions, n_samples):
    rho = scipy.stats.loguniform.rvs(a, b, size=(n_samples))
    # rho = np.random.uniform(0, 1, n_samples) ** (1 / d)
    r = np.abs(np.random.normal(size=(n_dimensions, n_samples)))
    c = np.sqrt(np.sum(r ** 2, axis=0))
    return rho * r / c


# Below is a 2-dim example

# Hypersphere
x = hypersphere(n_dimensions=2, n_samples=1000)

plt.figure()
plt.scatter(x[0]*500, x[1]*500)
plt.axis("equal")
plt.show()

# Log hypersphere

p_a = -3
p_b = 3

# a is min and b is max
y = log_hypersphere(10 ** p_a, 10 ** p_b, n_dimensions=8, n_samples=1000)

# Verify uniform distribution of different radia
print("Distribution of radia for log hypershere")
r = np.sqrt(np.sum(y ** 2, axis=0))
for i in range(p_a, p_b):
    print(f"Radius: {10**i} to {10**(i+1)}")
    print(len(r[(r < 10 ** (i + 1)) & (r > 10 ** i)]))

plt.figure()
plt.scatter(y[0], y[1])
plt.axis("equal")
plt.show()
# %%

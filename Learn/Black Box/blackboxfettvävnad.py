# Göra blackbox för fettvävnad
# Data från expdata_GLUCOSE_dr

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Definiera den nya funktionen??
def func(x, a, b, c, d):
    return a * np.exp(x * -b) + c + d

# Glukosupptag % av max 
ydata = np.array([0., 22.4816, 52.6589, 70.9468, 84.0467, 89.6239, 95.3122, 100])

# Isulinkocentration
xdata = np.array([1e-12, 1e-11, 3e-11, 1e-10, 3e-10, 1e-9, 1e-8, 1e-7])
plt.plot(xdata, ydata, 'b-', label='data')

# Fit for the parameters a, b, c of the function func
popt, pcov = curve_fit(func, xdata, ydata)
popt
A1 = np.array([ 2.554237065,  1.35190947,  0.47450618,  7.84894589])
plt.plot(xdata, func(xdata, *popt), 'r-',
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f' % tuple(popt))

# Constrain the optimization to the region of 0 <= a <= 3, 0 <= b <= 1 and 0 <= c <= 0.5:
popt, pcov = curve_fit(func, xdata, ydata, bounds=(-1000000, [1000000., 1000000., 1000000., 1000000.]))
popt
A2 = np.array([ 2.43708906,  1.45354748,  0.35015434,  8.88473444])
plt.semilogx(xdata, func(xdata, *popt), 'g--',
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f' % tuple(popt))

plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()

# Logaritmisk skala på x-axeln
#plt.semilogx(I, G)

# Plotta

#plt.xlabel("Insulin (M)", fontsize=12), plt.ylabel("Glukosupptag % av max", fontsize=12)
#plt.title("Glukosupptag")

#plt.show()
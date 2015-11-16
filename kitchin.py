from scipy.optimize import leastsq
import numpy as np

def loadValues(filename):
    'Loads a tab txt file with volumes and energies into numpy arrays.'
    values = np.genfromtxt(filename)

    vols_array = []
    energies_array = []
    for pair in values:
        vols_array.append(pair[0])
        energies_array.append(pair[1])

    vols = np.array(vols_array)
    energies = np.array(energies_array)

    return vols, energies

def Murnaghan(parameters, vol):
    'From Phys. Rev. B 28, 5480 (1983)'
    E0, B0, BP, V0 = parameters

    E = E0 + B0 * vol / BP * (((V0 / vol)**BP) / (BP - 1) + 1) - V0 * B0 / (BP - 1.0)

    return E

def objective(pars, y, x):
    #we will minimize this function
    err =  y - Murnaghan(pars, x)
    return err


x0 = [ -56.0, 0.54, 2.0, 16.5] #initial guess of parameters

vols, energies = loadValues('example.txt')

plsq = leastsq(objective, x0, args=(energies, vols))

print 'Fitted parameters = {0}'.format(plsq[0])

import matplotlib.pyplot as plt
plt.plot(vols,energies, 'ro')

#plot the fitted curve on top
x = np.linspace(min(vols), max(vols), 50)
y = Murnaghan(plsq[0], x)
plt.plot(x, y, 'k-')
plt.xlabel('Volume')
plt.ylabel('Energy')
plt.savefig('images/nonlinear-curve-fitting.png')

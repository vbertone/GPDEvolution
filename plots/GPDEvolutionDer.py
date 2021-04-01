import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings
from matplotlib.ticker import MaxNLocator
from matplotlib import ticker

# Loada data
data = np.loadtxt("GPDEvolutionDer.dat")

nd = 1000

x = np.reshape(data[:,1], (nd, nd))
mu = np.reshape(data[:,0], (nd, nd))
der = np.reshape(abs(data[:,2]), (nd, nd))

plt.title(r"\textbf{LO evolution from $\mu_0=1$ GeV, $\xi = 0.5$}", fontsize = 20)

plt.ylabel(r"\textbf{$\mu$ [GeV]}")
plt.xlabel(r"\textbf{$x$}")
plt.xlim(0, 1)
plt.ylim(1, 1000)
plt.yscale("log")
plt.vlines(0.5, 1, 1000, ls = "--", lw = 0.5)
plt.text(0.05, 500, r"$\displaystyle \left|\frac{dH_u(x,\xi,\mu)}{d\ln\mu^2}\right|$", fontsize = 18)

levels = np.logspace(-7, 1, 17)
img=plt.contourf(x, mu, der, levels = levels, locator = ticker.LogLocator(), cmap = "Spectral")
cbar = plt.colorbar(img, ticks = [0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10])

plt.savefig("GPDEvolutionDer.pdf")
plt.close()

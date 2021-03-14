import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

# Loada data
data = np.loadtxt("GPDMomentsPion.dat")

plt.title(r"\textbf{$\pi$ GPD, LO evolution from $\mu_0=1$ GeV to $\mu = 10$ GeV}", fontsize = 18)
#plt.text(0.0012, 0.28, r"\textbf{MMHT2014lo68cl}", fontsize = 16)

plt.ylabel(r"\textbf{Moments}", fontsize = 18)
plt.xlabel(r"$\xi$")
plt.xlim([-0.1, 1.1])
plt.ylim([0.4, 1.6])
plt.plot(data[:,0], data[:,1], "ro", label = r"\textbf{$\displaystyle \int_0^1 dx\,H_u^{(-)}(x,\xi,\mu)$}")
plt.plot(data[:,0], data[:,2], "bo", label = r"\textbf{$\displaystyle \int_0^1 dx\,x\left[\sum_q H_q^{(+)}(x,\xi,\mu)+H_g(x,\xi,\mu)\right]$}")

xi = np.linspace(0, 1, 100)
plt.plot(xi, (1+xi**2)/2, 'b--', lw = 1.5, label = r"$\displaystyle f(\xi) = \frac{1+\xi^2}2$")

plt.legend(fontsize = 14)
plt.savefig("GPDMomentsPion.pdf")
plt.close()

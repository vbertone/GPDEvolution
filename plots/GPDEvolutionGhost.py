import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

# Loada data
data = np.loadtxt("GPDEvolutionGhost.dat")

plt.title(r"\textbf{LO evolution from $\mu_0=1$ GeV, $\xi=0.1$}", fontsize = 20)

plt.ylabel(r"$H_u(x,\xi,\mu)$", fontsize = 18)
plt.xlabel(r"\textbf{$x$}")
plt.xlim([0, 1])
plt.ylim([-0.000000025, 0.000000025])
plt.plot(data[:,0], data[:,1] / data[:,0], "k--", lw = 1.5, label = r"\textbf{$\mu=1$ GeV}")
plt.plot(data[:,0], data[:,2] / data[:,0], "r-", label = r"\textbf{$\mu=\sqrt{10}$ GeV}")
plt.legend(fontsize = 18)

plt.savefig("GPDEvolutionGhost.pdf")
plt.close()

import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

# Loada data
data = np.loadtxt("ERBLEvolution.dat")

f, (ax1, ax2) = plt.subplots(2, 1, sharex = "all", gridspec_kw = dict(width_ratios = [1], height_ratios = [3, 1]))
plt.subplots_adjust(wspace = 0, hspace = 0)

ax1.set_title(r"\textbf{LO ERBL evolution from $\mu_0=1$ GeV}", fontsize = 20)

ax1.set_ylabel(r"$xF_{4}(x,\mu)$", fontsize = 18)
ax1.set_xlim([0.01, 1])
ax1.set_xscale("log")
ax1.set_ylim([-1.2, 1.2])
ax1.plot(data[:,0], data[:,1], label = r"\textbf{$\mu=1$ GeV}",   color = "black",  ls = "--", lw = 1.5)
ax1.plot(data[:,0], data[:,3], label = r"\textbf{$\mu=2$ GeV}",   color = "red",   )
ax1.plot(data[:,0], data[:,5], label = r"\textbf{$\mu=5$ GeV}",   color = "blue",  )
ax1.plot(data[:,0], data[:,7], label = r"\textbf{$\mu=10$ GeV}",  color = "orange",)
ax1.plot(data[:,0], data[:,9], label = r"\textbf{$\mu=100$ GeV}", color = "green", )

ax1.legend(fontsize = 18)

ax2.set_xlabel(r"$x$")
ax2.set_ylabel(r"\textbf{Ratio to analytic}", fontsize = 14)
ax2.set_xlim([0.01, 1])
ax2.set_xscale("log")
ax2.set_ylim([0.99, 1.01])
ax2.plot(data[:,0], data[:,1] / data[:,2],  color = "black", ls = "--", lw = 1.5)
ax2.plot(data[:,0], data[:,3] / data[:,4],  color = "red",   )
ax2.plot(data[:,0], data[:,5] / data[:,6],  color = "blue",  )
ax2.plot(data[:,0], data[:,7] / data[:,8],  color = "orange",)
ax2.plot(data[:,0], data[:,9] / data[:,10], color = "green", )

plt.savefig("ERBLEvolution.pdf")
plt.close()

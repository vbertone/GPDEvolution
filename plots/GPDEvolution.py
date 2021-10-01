import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

# Loada data
data = np.loadtxt("GPDEvolution.dat")

f, (ax1, ax2) = plt.subplots(2, 1, sharex = "all", gridspec_kw = dict(width_ratios = [1], height_ratios = [3, 1]))
plt.subplots_adjust(wspace = 0, hspace = 0)

ax1.set_title(r"\textbf{LO evolution from $\mu_0=1$ GeV to $\mu = 10$ GeV}", fontsize = 20)
ax1.text(0.0012, 0.5, r"\textbf{MMHT2014lo68cl}", fontsize = 16)

ax1.set_ylabel(r"$xF_u^{-}(x,\xi,\mu)$", fontsize = 18)
ax1.set_xlim([0.001, 1])
ax1.set_xscale("log")
ax1.set_ylim([-0.1, 1])
ax1.plot(data[:,0], data[:,2], label = r"\textbf{$\xi=0$}",    color = "red")
#ax1.plot(data[:,0], data[:,3], label = r"\textbf{$\xi=0.02$}", color = "blue")
ax1.plot(data[:,0], data[:,4], label = r"\textbf{$\xi=0.05$}", color = "blue")
#ax1.plot(data[:,0], data[:,5], label = r"\textbf{$\xi=0.2$}",  color = "green")
ax1.plot(data[:,0], data[:,6], label = r"\textbf{$\xi=0.5$}",  color = "orange")
#ax1.plot(data[:,0], data[:,7], label = r"\textbf{$\xi=0.9$}",  color = "cyan")
ax1.plot(data[:,0], data[:,8], label = r"\textbf{$\xi=1$}",    color = "green")
ax1.legend(fontsize = 18)

ax2.set_xlabel(r"$x$")
ax2.set_ylabel(r"\textbf{Ratio to DGLAP}", fontsize = 14)
ax2.set_xlim([0.001, 1])
ax2.set_xscale("log")
ax2.set_ylim([0, 3])
ax2.plot(data[:,0], data[:,1] / data[:,1], color = "red", ls = "--", lw = 1.5)
ax2.plot(data[:,0], data[:,2] / data[:,1], color = "red")
#ax2.plot(data[:,0], data[:,3] / data[:,1], color = "blue")
ax2.plot(data[:,0], data[:,4] / data[:,1], color = "blue")
#ax2.plot(data[:,0], data[:,5] / data[:,1], color = "green")
ax2.plot(data[:,0], data[:,6] / data[:,1], color = "orange")
#ax2.plot(data[:,0], data[:,7] / data[:,1], color = "cyan")
ax2.plot(data[:,0], data[:,8] / data[:,1], color = "green")

plt.savefig("GPDEvolution.pdf")
plt.close()

##############################################################
# Loada data
data = np.loadtxt("GPDEvolutionGluon.dat")

f, (ax1, ax2) = plt.subplots(2, 1, sharex = "all", gridspec_kw = dict(width_ratios = [1], height_ratios = [3, 1]))
plt.subplots_adjust(wspace = 0, hspace = 0)

ax1.set_title(r"\textbf{LO evolution from $\mu_0=1$ GeV to $\mu = 10$ GeV}", fontsize = 20)
ax1.text(0.0012, 0.5, r"\textbf{MMHT2014lo68cl}", fontsize = 16)

ax1.set_ylabel(r"$xF_g(x,\xi,\mu)$", fontsize = 18)
ax1.set_xlim([0.001, 1])
ax1.set_xscale("log")
ax1.set_yscale("log")
#ax1.set_ylim([-0.1, 1])
ax1.plot(data[:,0], data[:,2],  label = r"\textbf{$\xi=0$}",    color = "red")
ax1.plot(data[:,0], data[:,4],  label = r"\textbf{$\xi=0.05$}", color = "blue")
ax1.plot(data[:,0], data[:,8],  label = r"\textbf{$\xi=0.5$}",  color = "orange")
ax1.plot(data[:,0], data[:,11], label = r"\textbf{$\xi=1$}",    color = "green")
ax1.legend(fontsize = 18)

ax2.set_xlabel(r"$x$")
ax2.set_ylabel(r"\textbf{Ratio to DGLAP}", fontsize = 14)
ax2.set_xlim([0.001, 1])
ax2.set_xscale("log")
ax2.set_ylim([0, 3])
ax2.plot(data[:,0], data[:,1]  / data[:,1], color = "red", ls = "--", lw = 1.5)
ax2.plot(data[:,0], data[:,2]  / data[:,1], color = "red")
ax2.plot(data[:,0], data[:,4]  / data[:,1], color = "blue")
ax2.plot(data[:,0], data[:,8]  / data[:,1], color = "orange")
ax2.plot(data[:,0], data[:,11] / data[:,1], color = "green")

plt.savefig("GPDEvolutionGluon.pdf")
plt.close()

##############################################################
# Loada data
data = np.loadtxt("GPDEvolutionPlus.dat")

f, (ax1, ax2) = plt.subplots(2, 1, sharex = "all", gridspec_kw = dict(width_ratios = [1], height_ratios = [3, 1]))
plt.subplots_adjust(wspace = 0, hspace = 0)

ax1.set_title(r"\textbf{LO evolution from $\mu_0=1$ GeV to $\mu = 10$ GeV}", fontsize = 20)
ax1.text(0.0012, 0.5, r"\textbf{MMHT2014lo68cl}", fontsize = 16)

ax1.set_ylabel(r"$xF_u^+(x,\xi,\mu)$", fontsize = 18)
ax1.set_xlim([0.001, 1])
ax1.set_xscale("log")
ax1.set_ylim([-0.2, 1.8])
ax1.plot(data[:,0], data[:,2],  label = r"\textbf{$\xi=0$}",    color = "red")
ax1.plot(data[:,0], data[:,4],  label = r"\textbf{$\xi=0.05$}", color = "blue")
ax1.plot(data[:,0], data[:,8],  label = r"\textbf{$\xi=0.5$}",  color = "orange")
ax1.plot(data[:,0], data[:,11], label = r"\textbf{$\xi=1$}",    color = "green")
ax1.legend(fontsize = 18)

ax2.set_xlabel(r"$x$")
ax2.set_ylabel(r"\textbf{Ratio to DGLAP}", fontsize = 14)
ax2.set_xlim([0.001, 1])
ax2.set_xscale("log")
ax2.set_ylim([0, 3])
ax2.plot(data[:,0], data[:,1]  / data[:,1], color = "red", ls = "--", lw = 1.5)
ax2.plot(data[:,0], data[:,2]  / data[:,1], color = "red")
ax2.plot(data[:,0], data[:,4]  / data[:,1], color = "blue")
ax2.plot(data[:,0], data[:,8]  / data[:,1], color = "orange")
ax2.plot(data[:,0], data[:,11] / data[:,1], color = "green")

plt.savefig("GPDEvolutionPlus.pdf")
plt.close()

##############################################################
# Loada data
data = np.loadtxt("GPDMoments.dat")

plt.title(r"\textbf{LO evolution from $\mu_0=1$ GeV to $\mu = 10$ GeV}", fontsize = 20)
plt.text(0.7, 3, r"\textbf{MMHT2014lo68cl}", fontsize = 16)

plt.ylabel(r"$\displaystyle \int_0^1 dx\,x^{2n} F_u^{-}(x,\xi,\mu)$", fontsize = 18)
plt.xlabel(r"$\xi$")
plt.xlim([-0.1, 1.1])
plt.ylim([0.04, 10])
plt.yscale("log")
plt.plot(data[:,0], data[:,1], "ro", label = r"\textbf{$n=0$}")
plt.plot(data[:,0], data[:,3], "bo", label = r"\textbf{$n=1$}")

# Fit with plynomials of the appropriate degree
p0 = np.poly1d(np.polyfit(data[:,0], data[:,1], 0))
p2 = np.poly1d(np.polyfit(data[:,0], data[:,3], 2))
print(np.polyfit(data[:,0], data[:,1], 0))
print(np.polyfit(data[:,0], data[:,3], 2))

xi = np.linspace(0, 1, 100)
plt.plot(xi, p0(xi), 'r--', lw = 1.5, label = r"\textbf{Fit with $f(\xi)=A_0$}")
plt.plot(xi, p2(xi), 'b--', lw = 1.5, label = r"\textbf{Fit with $f(\xi)=A_0+A_1\xi^2$}")

plt.legend(fontsize = 18, ncol = 2)
plt.savefig("GPDMoments.pdf")
plt.close()

##############################################################
# Loada data
data = np.loadtxt("GPDMoments.dat")

plt.title(r"\textbf{LO evolution from $\mu_0=1$ GeV to $\mu = 10$ GeV}", fontsize = 20)
plt.text(0.7, 0.4, r"\textbf{MMHT2014lo68cl}", fontsize = 16)

plt.ylabel(r"$\displaystyle \int_0^1 dx\,x^{2n+1} F_u^{+}(x,\xi,\mu)$", fontsize = 18)
plt.xlabel(r"$\xi$")
plt.xlim([-0.1, 1.1])
plt.ylim([0.015, 1])
plt.yscale("log")
plt.plot(data[:,0], data[:,2], "ro", label = r"\textbf{$n=0$}")
plt.plot(data[:,0], data[:,4], "bo", label = r"\textbf{$n=1$}")

# Fit with plynomials of the appropriate degree
p1 = np.poly1d(np.polyfit(data[:,0], data[:,2], 2))
p3 = np.poly1d(np.polyfit(data[:,0], data[:,4], 4))
print(np.polyfit(data[:,0], data[:,2], 2))
print(np.polyfit(data[:,0], data[:,4], 4))

xi = np.linspace(0, 1, 100)
plt.plot(xi, p1(xi), 'r--', lw = 1.5, label = r"\textbf{Fit with $f(\xi)=B_0 + B_1\xi^2$}")
plt.plot(xi, p3(xi), 'b--', lw = 1.5, label = r"\textbf{Fit with $f(\xi)=B_0 + B_1\xi^2 + B_2\xi^4$}")

plt.legend(fontsize = 18, ncol = 2)
plt.savefig("GPDMomentsPlus.pdf")
plt.close()













##############################################################
# Loada data
data = np.loadtxt("GPDMomentsGK.dat")

plt.title(r"\textbf{LO evolution from $\mu_0=\sqrt{2}$ GeV to $\mu = 10$ GeV}", fontsize = 20)
plt.text(0.7, 3, r"\textbf{GK16 model}", fontsize = 16)

plt.ylabel(r"$\displaystyle \int_0^1 dx\,x^{2n} F_u^{-}(x,\xi,\mu)$", fontsize = 18)
plt.xlabel(r"$\xi$")
plt.xlim([-0.1, 1.1])
plt.ylim([0.04, 10])
plt.yscale("log")
plt.plot(data[:,0], data[:,1], "ro", label = r"\textbf{$n=0$}")
plt.plot(data[:,0], data[:,3], "bo", label = r"\textbf{$n=1$}")

# Fit with plynomials of the appropriate degree
p0 = np.poly1d(np.polyfit(data[:,0], data[:,1], 0))
p2 = np.poly1d(np.polyfit(data[:,0], data[:,3], 2))
print(np.polyfit(data[:,0], data[:,1], 0))
print(np.polyfit(data[:,0], data[:,3], 2))

xi = np.linspace(0, 1, 100)
plt.plot(xi, p0(xi), 'r--', lw = 1.5, label = r"\textbf{Fit with $f(\xi)=A_0$}")
plt.plot(xi, p2(xi), 'b--', lw = 1.5, label = r"\textbf{Fit with $f(\xi)=A_0+A_1\xi^2$}")

plt.legend(fontsize = 18, ncol = 2)
plt.savefig("GPDMomentsGK.pdf")
plt.close()

##############################################################
# Loada data
data = np.loadtxt("GPDMomentsGK.dat")

plt.title(r"\textbf{LO evolution from $\mu_0=\sqrt{2}$ GeV to $\mu = 10$ GeV}", fontsize = 20)
plt.text(0.7, 0.4, r"\textbf{GK16 model}", fontsize = 16)

plt.ylabel(r"$\displaystyle \int_0^1 dx\,x^{2n+1} F_u^{+}(x,\xi,\mu)$", fontsize = 18)
plt.xlabel(r"$\xi$")
plt.xlim([-0.1, 1.1])
plt.ylim([0.015, 1])
plt.yscale("log")
plt.plot(data[:,0], data[:,2], "ro", label = r"\textbf{$n=0$}")
plt.plot(data[:,0], data[:,4], "bo", label = r"\textbf{$n=1$}")

# Fit with plynomials of the appropriate degree
p1 = np.poly1d(np.polyfit(data[:,0], data[:,2], 2))
p3 = np.poly1d(np.polyfit(data[:,0], data[:,4], 4))
print(np.polyfit(data[:,0], data[:,2], 2))
print(np.polyfit(data[:,0], data[:,4], 4))

xi = np.linspace(0, 1, 100)
plt.plot(xi, p1(xi), 'r--', lw = 1.5, label = r"\textbf{Fit with $f(\xi)=B_0 + B_1\xi^2$}")
plt.plot(xi, p3(xi), 'b--', lw = 1.5, label = r"\textbf{Fit with $f(\xi)=B_0 + B_1\xi^2 + B_2\xi^4$}")

plt.legend(fontsize = 18, ncol = 2)
plt.savefig("GPDMomentsGKPlus.pdf")
plt.close()


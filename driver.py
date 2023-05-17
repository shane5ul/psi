import numpy as np
import psi as prony

# Single mode Maxwell with g1 = 10.0, tau1 = 10.0, and Ge = 1.0
g   = np.array([10.0, 1.0])
tau = np.array([10.0])

# use PSI to compute j, lam
j, lam = prony.findJtProny(g, tau)

# use computed j, lam to interconvert back
g, tau = prony.findGtProny(j, lam)

# to plot
t = np.geomspace(1e-2, 1e2)
G = prony.evalGtProny(t, g, tau, isPlot=True)
J = prony.evalJtProny(t, j, lam, isPlot=True)

# using file: Prony series
# savefile=True prints to "retardSpect.dat" or "relaxSpect.dat" as appropriate
j, lam = prony.fileG2J('tests/g2.dat', savefile=False)
g, tau = prony.fileJ2G('tests/l4.dat', savefile=False)

#
#  Program: PSI (Prony Series Interpolation)
#  Date   : May 15, 2023
#
# Convention adopted
# ------------------
# g = g for liq, and [g, Ge] for solids, where Ge is the plateau modulus
# j = [j Je invEta0] for liq, and [j Je] for solids
#
#  Usage  :
# 
#        If you have arrays (g, tau) or (j, lam) then for interconversion use:
#        j, lam = findJtProny(g, tau)
#        g, tau = findGtProny(j, lam)
#    
#        If you have a file like dmodes.dat generated from pyReSpect-time or pyReSpect-Jt
#        fileJ2G(fname, savefile)
#        fileG2J(fname, savefile)
#        fname points to the file containing the modes;
#           the first line should have # Ge, Je, invEta0 info if nonzero (otherwise)
#              # Je = 1.0 [viscoelastic solid; relaxation spectrum]
               # Ge = 1.0 [viscoelastic solid; retardation spectrum]
               # Je, invEta0 = 7.3970e-05 9.9080e-05 [viscoelastic liq; retardation spectrum]
#        savefile = True saves the interconverted prony series to relaxSpect.dat or retardSpect.dat
#
#        To evaluate the Prony series use:
#        G = evalGtProny(t, g, tau, isPlot)
#        J = evalJtProny(t, j, lam, isPlot)
#


import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from scipy.signal import invres

#plt.style.use(['myjournal', 'seaborn-ticks'])

########### H E L P E R   F U N C T I O N S ###########

def evalJtProny(t, j, lam, isPlot=False):
    """
    This evaluates J(t) given the DRS j = [j1, ..., jnmodes, Je, invEta0]
    """
    
    nmodes = len(lam)
    nex    = len(j) - nmodes # 1 for solid, 2 for liq
    
    Je = j[nmodes]
    
    if nex > 1:
        invEta0 = j[nmodes+1]
    else:
        invEta0 = 0.
    
    # given coefficients, evaluate Jt
    Jt = Je + invEta0 * t
    for i in range(nmodes):
        Jt -= j[i] * np.exp(-t/lam[i])

    if isPlot:    
        plt.plot(t, Jt)

        for i in range(nmodes):
            plt.axvline(x=lam[i], c='gray', alpha=0.5)
        
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$J(t)$')
        plt.tight_layout()
        plt.show()    
        
    return Jt

def evalGtProny(t, g, tau, isPlot=False):
    """
    This evaluates G(t) given the DRS g = [g1, ..., g_nmodes, Ge]
    """

    G      = np.zeros(len(t))
    nmodes = len(tau)

    # solid or liquid
    if len(g) > nmodes:
        G += g[nmodes] # N = N_s + 1
        
    for i in range(nmodes):
        G += g[i] * np.exp(-t/tau[i])
    
    if isPlot:
        plt.plot(t, G)
        plt.ylim(max(np.min(G),np.amax(G)/1e6), None)

        for i in range(nmodes):
            plt.axvline(x=tau[i], c='gray', alpha=0.5)

        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$G(t)$')
        plt.tight_layout()
        plt.show() 

    return G

def readRelaxSpec(fname):
    """
    read input file containing relaxation spectrum
    Filename fname contains (g, tau):
    header contains Ge on first line # Ge = xxxx   
    """

    # Read header
    try:
        with open(fname) as f:
            first_line = f.readline().split()
    except OSError:
        print('Problem reading file ' + fname)
        quit()

    if 'Je' in first_line:
        print("Input file " + fname + " contains Je in header: use readRetardSpec?")
        quit()


    # can read body with two or three columns
    data = np.loadtxt(fname, ndmin=2)

    g    = data[:,0]; 
    tau  = data[:,1] 

    if np.isscalar(tau):
        tau = np.array([tau])
        g   = np.array([g])
        
    nmode = len(tau)

    if 'Ge' in first_line:
        Ge      = float(first_line[-1]);
        g       = np.append(g, Ge)
    
    return g, tau
    
def readRetardSpec(fname):
    """
    Filename fname contains (j, lambda):
    header contains Je [and optionally invEta0] on first line
    """

    # Read Je and invEta0
    try:
        with open(fname) as f:
            first_line = f.readline().split()
    except OSError:
        print('Problem reading file ' + fname)
        quit()

    if 'Ge' in first_line:
        print("Input file " + fname + " contains Ge in header: use readRelaxSpec?")
        quit()

    # can read body with two or three columns (ndmin protects when only one mode)
    data = np.loadtxt(fname, ndmin=2)
    j    = data[:,0]; 
    lam  = data[:,1] 

    if np.isscalar(lam):
        lam = np.array([lam])
        j   = np.array([j])

    nmode     = len(lam)

    if 'invEta0' in first_line:
        isLiq = True
        invEta0 = float(first_line[-1]); 
        Je      = float(first_line[-2]);
        j       = np.append(j, [Je, invEta0])
    else:
        isLiq = False
        Je      = float(first_line[-1]);
        j       = np.append(j, Je)


    # Je ~ sum(j), but writing errors leading to -ve
    if (1 - sum(j[:nmode])/Je) < 1e-4:
        j[nmode] = sum(j[:nmode])

    return j, lam

#########  P S I   C O R E  F U N C T I O N S #########

def findBeta(g, tau, isPlot=False):
    """
    find the roots -beta of G(s); return positive rates +beta

    4/13/2023: using new method of first finding the numerator for additional stability
    and absence of poles

    no longer directly using interpolation property: since roots seems to work fine!
    + accounts naturally when polynomial is one degree smaller than expected
    """    
    
    nmodes = len(tau)
    alpha = 1.0 / tau
    # if solid then alpha0 = 0; append at end because that is where Ge sits
    if len(g) > nmodes: 
        alpha = np.append(alpha, 0.0)

        
    # this builds Gs into polynomials with coefficients num/den
    num, _  = invres(g, -alpha, [], tol=1e-12)  ##########################
    negBeta = np.roots(num)                     #

    # older approach using interpolation property: but leads to the same roots at np.roots
    # negBeta  = np.zeros(len(alpha) - 1) # equal to #taus (solid) or 1 less (liq)
    #
    # numPoly = np.poly1d(num)    # polynomial function
    # # find the roots using bracketing: finding exactly the same roots!
    # for i in range(len(beta)):
    #     beta[i] = root_scalar(numPoly, bracket = [-alpha[i], -alpha[i+1]],
    #                         method='brentq').root     
    #     print('{:3d}\t{:.4e} {:.4e} {:.4e}\t{:.4e}'.format(i, -alpha[i], beta[i], -alpha[i+1], numPoly(beta[i])))

    beta    = -negBeta

    if isPlot:
        plt.plot(-alpha, np.zeros(len(alpha)), 'o', alpha=0.3)
        plt.plot(beta, np.zeros(len(beta)), 'x')
        plt.xscale('log')
        plt.show()

    return beta

def findAlpha(j, lam, isPlot=False):
    """
    find the roots -alpha of J(s); return positive rates +alpha

    use new method of first finding the numerator for additional stability and absence of poles
    
    If Je = sum(j), the numerator drops one order: J(0) = 0 which is unphysical
    
    """
    nmodes = len(lam)
    beta   = 1.0/lam


    if len(j) - nmodes == 1: # solid    
        r = np.append(j[nmodes], -j[:nmodes])
        p = np.append(0., -beta)
        # alpha = np.zeros(len(beta))
    else:
        r = np.append([j[nmodes], j[nmodes+1]], -j[:nmodes])
        p = np.append([0., 0.], -beta)
        # alpha = np.zeros(len(beta)+1)
        
    num, den = invres(r, p, [], tol=1e-12)
    negAlpha = np.roots(num)
    alpha    = -negAlpha 

    if isPlot:
        plt.plot(alpha, np.zeros(len(alpha)), 'o', alpha=0.3, label=r'$\alpha$')
        plt.plot(beta, np.zeros(len(beta)), 'x', label=r'$\beta$')
        plt.xscale('log')
        plt.legend()
        plt.show()
        
    return alpha

def findJtProny(g, tau):
    """
    takes in (g, tau) returns (j, lam)
    """
    
    nmodes = len(tau)
    
    # sort the modes if nmodes > 1
    if nmodes > 1:
        idx    = np.argsort(tau)
        tau    = tau[idx]
        g[:nmodes] = g[:nmodes][idx]
    
    beta   = findBeta(g, tau)
    alpha  = 1.0 / tau    
    j      = np.zeros(len(beta))

    # solid or liquid
    if len(g) > nmodes: 
        Je    = 1/g[nmodes]
        alpha = np.append(alpha, 0.0)
        isLiq = False    
    else:
        invEta0 = 1.0/np.sum(g * tau)
        Je      = invEta0**2 * np.sum(g * tau**2)
        isLiq   = True
    
    # find the coefficients of the Jt Prony series
    for i in range(len(j)):
        j[i] = 1.0 / np.sum(g*beta[i]**2 / (alpha - beta[i])**2)

    # pack appropriately
    jPlus = np.append(j, Je)    
    if isLiq:
        jPlus = np.append(jPlus, invEta0)
        
    return jPlus, 1.0/beta 

def findGtProny(j, lam):
    """
    takes (j, lam) and returns (g, tau)
    """
    nmodes = len(lam)

    # sort the modes by lambda
    idx    = np.argsort(lam)
    lam    = lam[idx]
    j[:nmodes] = j[:nmodes][idx]

    alpha  = findAlpha(j, lam)
    beta   = 1.0/lam
    g      = np.zeros(len(alpha))

    # solid or liquid
    if len(j) == nmodes+1:
        Je      = j[-1]
        invEta0 = 0.
        isLiq   = False
        Ge      = 1/Je
    else:
        Je      = j[-2]
        invEta0 = j[-1]
        isLiq = True

    # find the coefficients of the Gt Prony series
    for i in range(len(g)):
        tmp  = np.sum(j[:nmodes]*alpha[i]**2 / (beta - alpha[i])**2)
        g[i] = 1.0 /(-Je + 2*invEta0/alpha[i] + tmp)
        
        
    # pack appropriately
    if not isLiq:
        g = np.append(g, Ge)
        
    return g, 1.0/alpha

def fileJ2G(fname, savefile=True):
    """
    1. read file
    2. run interconv calc
    3. save result
    """

    j, lam = readRetardSpec(fname)
    g, tau = findGtProny(j, lam)
    
    N = len(tau)
        
    if savefile:
        if len(g) > N:
            np.savetxt('relaxSpect.dat', np.c_[g[:N], tau], fmt='%e', 
                        header='Ge = {0:0.4e}'.format(g[N]))        
        else:
            np.savetxt('relaxSpect.dat', np.c_[g, tau])
        
    return g, tau

def fileG2J(fname, savefile=True):
    """
    1. read file
    2. run interconv calc
    3. save result
    """
    g, tau = readRelaxSpec(fname)

    j, lam = findJtProny(g, tau)
    N = len(lam)
    if savefile:
        if len(g) > len(tau):
            np.savetxt('retardSpect.dat', np.c_[j[:N], lam], fmt='%e', 
                        header='Je = {0:0.4e}'.format(j[N]))        
        else:
            np.savetxt('retardSpect.dat', np.c_[j[:N], lam], fmt='%e', 
                        header='Je, invEta0 = {0:0.4e} {1:0.4e}'.format(j[N], j[N+1]))        
        
    return j, lam

#    
# Main Driver: This part is not run when imported as a module
# Otherwise you can directly include commands below the following "if" statement.
#
if __name__ == '__main__':

    # G -> J
    j, lam = fileG2J('tests/g2.dat', savefile=True)

    t = np.geomspace(1e-3, 1e3)
    J = evalJtProny(t, j, lam, isPlot=True)

    # J -> G
    g, tau = fileJ2G('tests/l1.dat', savefile=True)
    G = evalGtProny(t, g, tau, isPlot=True)
    
    # if you know DRS then don't need to access fileG2J and fileJ2G
    j, lam = findJtProny(g, tau)
    g, tau = findGtProny(j, lam)

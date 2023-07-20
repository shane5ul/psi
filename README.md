# PSI 
PSI (**Prony Series Interconversion**) is a computer program for interconversion of relaxation modulus $G(t)$ and creep compliance $J(t)$. It is based on the method of Loy et al. (2015) with a few modifications:

+ Shanbhag, "A computer program for interconversion of creep compliance and relaxation modulus", Journal of Rheology 67, 965–975 **2023** [doi: 10.1122/8.0000695]
+ Loy, de Hoog, and Anderssen, "Interconversion of Prony series for relaxation and creep", Journal of Rheology 59, 1261 **2015** [doi: 10.1122/1.4929398]

## Prony Series

The modulus and compliance are represented using a Prony series. The convention used to pack the coefficients and time constants of the Prony series in PSI is also provided.

**Solids**

$$G(t) = G_e + \sum\limits_{k=1}^{N} g_k e^{-t/\tau_k}$$

$$J(t) = J_e - \sum\limits_{k=1}^{N} j_{k} e^{-t/\lambda_{k}}$$

`g = ` $[g_1, g_2, \cdots g_N, G_e]$ and `tau = ` $[\tau_1, \tau_2, \cdots, \tau_N]$

`j = ` $[j_1, j_2, \cdots, j_N, J_e]$ and `lam = ` $[\lambda_1, \lambda_2, \cdots, \lambda_N]$

**Liquids**

$$G(t) = \sum\limits_{k=1}^{N} g_k e^{-t/\tau_k}$$

$$J(t) = J_e + \eta_0^{-1} t - \sum\limits_{k=1}^{N-1} j_{k} e^{-t/\lambda_{k}}$$

`g = ` $[g_1, g_2, \cdots g_N]$ and `tau = ` $[\tau_1, \tau_2, \cdots, \tau_N]$

`j = ` $[j_1, j_2, \cdots, j_{N-1}, J_e, \eta_{0}^{-1}]$ and `lam = ` $[\lambda_1, \lambda_2, \cdots, \lambda_{N-1}]$


## Code

The program `psi.py` contains all the necessary functions.

###  Usage  :

`psi.py` can be imported as a module to access the following functions.

Alternatively, instructions for interconversion can be directly included in `psi.py`. The program can be run using
`python3 psi.py`
 
If you have arrays `(g, tau)` or `(j, lam)` then for interconversion then use

`j, lam = findJtProny(g, tau)`

`g, tau = findGtProny(j, lam)`

If you have a file like the output `dmodes.dat` generated from `pyReSpect-time` or `pyJt` 

`j, lam = fileG2J(fname)`

`g, tau = fileJ2G(fname)`

where `fname` points to the file containing the modes. The first line should have # Ge, Je, invEta0 info if nonzero. See included examples.

`# Je = 1.0`

`# Ge = 1.0`

`# Je, invEta0 = 7.3970e-05 9.9080e-05`

Some example files are included in the tests/ folder.

If the Prony series $G(t)$ or $J(t)$ can be evaluated using,

`G = evalGtProny(t, g, tau, isPlot=True)`

`J = evalJtProny(t, j, lam, isPlot=True)`

## Pre-requisites

This code has been tested on a Linux box. The numbers in parenthesis show the version it has been tested on. 

python3 (3.6.9)
numpy (1.19.5)
scipy (1.5.4)
matplotlib (3.3.4)

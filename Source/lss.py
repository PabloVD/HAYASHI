""""
Large scale structure routines
Pablo Villanueva Domingo
Started: October 2018
Last update: September 2022
"""

from Source.cosmo import *
from colossus.lss import mass_function

# log cutoff in integrals
cutoff = 10.
# num of k values in k integrals, set 1000 for dndlnM plot
num_k = 1000

#--- Standard large scale structure ---#

# Radius conversion from a given mass
def MtoR(M):
    return (4.*np.pi/3.*Omega_m*rho_c_Mpc)**(-1./3.)*M**(1./3.)

# Top-hat window
def WindowTH(x):
    return 3*(np.sin(x)-x*np.cos(x))/x**3

# Logarithmic derivative of the top-hat window, dW/dlnR
def DerWindowTH(x):
    return (9*x*np.cos(x)+3*(x**2 - 3)*np.sin(x))/x**3

# Press-Schechter first crossing distribution
def f_ps(x):
    return np.sqrt(2.*x/np.pi)*np.exp(-x/2.)

# Sheth-Tormen first crossing distribution
def f_st(x):
    p = 0.3
    q = 0.707
    A = 0.3222
    return A*np.sqrt(2*q*x/np.pi)*( 1 + (q*x)**(-p) )*np.exp(-q*x/2)

# Halo mass function in LCDM cosmology, dndlnM, units of (h/Mpc)^3, using the Sheth-Tormen parameterization
def dndlnM_CDM(M,z):
    return mass_function.massFunction(M, z, q_out = "dndlnM", model = "sheth99")

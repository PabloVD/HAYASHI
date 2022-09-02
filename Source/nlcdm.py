""""
Cosmological functions
Pablo Villanueva Domingo
Started: October 2018
Last update: April 2021
"""

from Source.functions import *


#--- POWER SPECTRUM ---#

# Isocurvature power spectrum for PBHs at z=0, (Mpc/h)^3
def IsoPS(k,fpbh,Mpbh):
    onepluszeq = 2.4e4*Omega_m*hlittle**2
    keq = 0.073*Omega_m*hlittle    # h/Mpc
    if k > keq:
        return (9./4.)*onepluszeq**2.*fpbh*Mpbh*hlittle/((Omega_m-Omega_b)*rho_c_Mpc)
    else:
        return 0.

# Total power spectrum for PBH universes at z=0, (Mpc/h)^3
def PkPBH(k,fpbh,Mpbh):
    IsoPSvec = np.vectorize(IsoPS)
    return cosmo.matterPowerSpectrum(k) + IsoPSvec(k,fpbh,Mpbh)

# Radius conversion from a given mass
def MtoR(M):
    return (4.*np.pi/3.*Omega_m*rho_c_Mpc)**(-1./3.)*M**(1./3.)

# Top-hat window
def WindowTH(x):
    return 3*(np.sin(x)-x*np.cos(x))/x**3

# Logarithmic derivative of the top-hat window, dW/dlnR
def DerWindowTH(x):
    return (9*x*np.cos(x)+3*(x**2 - 3)*np.sin(x))/x**3

# Variance of the linear field in PBH scenarios
def Sigma2(M,fpbh,Mpbh):
    R = MtoR(M)
    logk = np.linspace(-cutoff,cutoff,num=num_k)
    return integrate.simps((1./(2.*np.pi**2.))*np.exp(3.*logk)*PkPBH(np.exp(logk),fpbh,Mpbh)*WindowTH(np.exp(logk)*R)**2, logk )
    #return integrate.quad(lambda logk: (1./(2.*np.pi**2.))*np.exp(3.*logk)*PkPBH(np.exp(logk),fpbh,Mpbh)*WindowTH(np.exp(logk)*R)**2, -cutoff, cutoff)[0]

# dsigma^2/dlnM
def DerSigma2(M,fpbh,Mpbh):
    R = MtoR(M)
    logk = np.linspace(-cutoff,cutoff,num=num_k)
    return integrate.simps((2./3.)*(1/(2.*np.pi**2.))*np.exp(3.*logk)*PkPBH(np.exp(logk),fpbh,Mpbh)*WindowTH(np.exp(logk)*R)*DerWindowTH(np.exp(logk)*R), logk )
    #return integrate.quad(lambda logk: (2./3.)*(1/(2.*np.pi**2))*np.exp(3.*logk)*PkPBH(np.exp(logk),fpbh,Mpbh)*WindowTH(np.exp(logk)*R)*DerWindowTH(np.exp(logk)*R), -cutoff, cutoff)[0]

#--- HALOS ---#

# Press-Schechter first crossing distribution
def f_ps(x):
    return np.sqrt(2.*x/np.pi)*np.exp(-x/2.)

# Sheth-Tormen first crossing distribution
def f_st(x):
    p = 0.3
    q = 0.707
    A = 0.3222
    return A*np.sqrt(2*q*x/np.pi)*( 1 + (q*x)**(-p) )*np.exp(-q*x/2)

# Halo mass function for PBH scenarios, dndlnM, units of (h/Mpc)^3
def dndlnM_PBH(M,z,fpbh,Mpbh):
    s2 = sigvec(M,fpbh,Mpbh)
    deltac = 1.68647/cosmo.growthFactor(z)
    nu = (deltac**2/s2)
    dlogsigdlogm = np.abs(dersigvec(M,fpbh,Mpbh)/(2.*s2))
    return Omega_m*rho_c_Mpc*f_st(nu)/M*dlogsigdlogm

# Halo mass function in LCDM cosmology, dndlnM, units of (h/Mpc)^3
def dndlnM_CDM(M,z):
    return mass_function.massFunction(M,z, q_out = "dndlnM", model = "sheth99")

# Vectorize some functions to allow them to take arrays as inputs
sigvec = np.vectorize(Sigma2)
dersigvec = np.vectorize(DerSigma2)
dndlnM_PBHvec = np.vectorize(dndlnM_PBH)
#OptDepth21 = np.vectorize(OptDepth21)

""""
Non-LambdaCDM large scale structure (warm dark matter and primordial black holes)
Pablo Villanueva Domingo
Started: October 2018
Last update: September 2022
"""

from Source.lss import *


#--- Primordial Black Holes ---#

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
    return cosmomodel.matterPowerSpectrum(k) + IsoPSvec(k,fpbh,Mpbh)

# Variance of the linear field in PBH scenarios
def Sigma2_PBH(M,fpbh,Mpbh):
    R = MtoR(M)
    logk = np.linspace(-cutoff,cutoff,num=num_k)
    return integrate.simps((1./(2.*np.pi**2.))*np.exp(3.*logk)*PkPBH(np.exp(logk),fpbh,Mpbh)*WindowTH(np.exp(logk)*R)**2, logk )

# dsigma^2/dlnM in PBH scenarios
def DerSigma2_PBH(M,fpbh,Mpbh):
    R = MtoR(M)
    logk = np.linspace(-cutoff,cutoff,num=num_k)
    return integrate.simps((2./3.)*(1/(2.*np.pi**2.))*np.exp(3.*logk)*PkPBH(np.exp(logk),fpbh,Mpbh)*WindowTH(np.exp(logk)*R)*DerWindowTH(np.exp(logk)*R), logk )

# Halo mass function for PBH scenarios, dndlnM, units of (h/Mpc)^3
def dndlnM_PBH(M,z,fpbh,Mpbh):
    s2 = Sigma2_PBH(M,fpbh,Mpbh)
    deltac = 1.68647/cosmomodel.growthFactor(z)
    nu = (deltac**2/s2)
    dlogsigdlogm = np.abs(DerSigma2_PBH(M,fpbh,Mpbh)/(2.*s2))
    return Omega_m*rho_c_Mpc*f_st(nu)/M*dlogsigdlogm



#--- Warm Dark Matter ---#

# WDM transfer function
def Transfer2(k,m):
    nu = 1.12
    alpha = 0.049*m**(-1.11)
    return (1.+(alpha*k)**(2.*nu))**(-10./nu)

# WDM Power spectrum
def ps_wdm(k,m):
    return Transfer2(k,m)*cosmomodel.matterPowerSpectrum(k)

# Half-mode length
def lamb_hmm(m):
    nu = 1.12
    alpha = 0.049*m**(-1.11)
    return 2.*np.pi*alpha*(2.**(nu/5.) - 1.)**(-1./2./nu)

# Half-mode mass
def M_hmm(m):
    return (4.*np.pi/3.*Omega_m*rho_c_Mpc)*(lamb_hmm(m)/2.)**3.

# Variance of the linear field in WDM scenarios
def Sigma2_WDM(M,mwdm):
    R = MtoR(M)
    logk = np.linspace(-cutoff,cutoff,num=num_k)
    return integrate.simps((1./(2.*np.pi**2.))*np.exp(3.*logk)*ps_wdm(np.exp(logk),mwdm)*WindowTH(np.exp(logk)*R)**2, logk )
    #return integrate.quad(lambda logk: (1./(2.*np.pi**2.))*np.exp(3.*logk)*PkPBH(np.exp(logk),fpbh,Mpbh)*WindowTH(np.exp(logk)*R)**2, -cutoff, cutoff)[0]

# dsigma^2/dlnM in WDM scenarios
def DerSigma2_WDM(M,mwdm):   # dsigma^2/dlnM
    R = MtoR(M)
    logk = np.linspace(-cutoff,cutoff,num=num_k)
    return integrate.simps((2./3.)*(1/(2.*np.pi**2.))*np.exp(3.*logk)*ps_wdm(np.exp(logk),mwdm)*WindowTH(np.exp(logk)*R)*DerWindowTH(np.exp(logk)*R), logk )
    #return integrate.quad(lambda logk: (2./3.)*(1/(2.*np.pi**2))*np.exp(3.*logk)*PkPBH(np.exp(logk),fpbh,Mpbh)*WindowTH(np.exp(logk)*R)*DerWindowTH(np.exp(logk)*R), -cutoff, cutoff)[0]

def dndlnM_WDM(M,z,mwdm):  # dndlnM units of (h/Mpc)^3
    s2 = Sigma2_WDM(M,mwdm)
    deltac = 1.68647/cosmomodel.growthFactor(z)
    nu = (deltac**2/s2)
    dlogsigdlogm = np.abs(DerSigma2_WDM(M,mwdm)/(2.*s2))
    return Omega_m*rho_c_Mpc*f_st(nu)/M*dlogsigdlogm*(1. + 2.*M_hmm(mwdm)/M)**(-0.6)

# Vectorize some functions to allow them to take arrays as inputs
Sigma2_PBH = np.vectorize(Sigma2_PBH)
DerSigma2_PBH = np.vectorize(DerSigma2_PBH)
Sigma2_WDM = np.vectorize(Sigma2_WDM)
DerSigma2_WDM = np.vectorize(DerSigma2_WDM)

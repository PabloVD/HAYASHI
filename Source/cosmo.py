""""
Cosmological functions
Pablo Villanueva Domingo
Started: October 2018
Last update: April 2021
"""


from scipy import integrate
from Source.constants import *

# Employ the Colossus package for some cosmological computations:
# https://bdiemer.bitbucket.io/colossus/tutorials.html
from colossus.cosmology import cosmology
from colossus.lss import mass_function
from colossus.halo.concentration import concentration



cosmomodel = cosmology.setCosmology('planck18')
rho_c_Mpc = cosmomodel.rho_c(0)*1e9   # in M_sun h^2 Mpc^-3

#--- COSMO ---#

# Hubble function, s^-1
def Hz(z):
    return H0*np.sqrt(Omega_m*(1.+z)**3. + Omega_l)

# dt/dz, in seconds
def dtdz(z):
    return -1./Hz(z)/(1.+z)

# dr/dz, in cm
def drdz(z):
    return -c*dtdz(z)*(1.+z)

# Adiabatic temperature, in K
def Tk_ad(z):
    return Tcmb0*(1.+zdec)*((1.+z)/(1.+zdec))**2.

# Jeans mass, Msun/h
def MJeans(z,T): # \simeq 3.58e5*(Tk_ad(z)/(1.+z))**(3./2.) Msun/hlittle ?
    lambdaJ = (5.*np.pi*T*kelvintoerg/(3.*G_N*rho_m_z(z)*mu*m_p))**(1./2.)
    #return 4.*np.pi*rho_m_z(z)/3.*(lambdaJ)**(3.)/(Msun/hlittle)   # \simeq 4.96e5*(Tk/(1.+z))**(3./2.) Msun
    return 4.*np.pi*rho_m_z(z)/3.*(lambdaJ/2.)**(3.)/(Msun/hlittle)   # \simeq 6.21e4*(Tk/(1.+z))**(3./2.) Msun
    #return 4.*np.pi*rho_m_z(z)/3.*(5.*3.*T*kelvintoerg/(4.*np.pi*G_N*rho_m_z(z)*mu*m_p))**(3./2.)/(Msun/hlittle)

# Matter density, g cm^-3
def rho_m_z(z):
    return Omega_m*rho_c*(1.+z)**3.

# Critical density as a function of z, g cm^-3
def rho_c_z(z):
    return rho_c*( Omega_m*(1.+z)**3. + Omega_l )

#--- HALO VIRIALITIES ---#

# Omega_m*(1.+z)**3./(Omega_m*(1.+z)**3. + Omega_l)
def Omz(z):
    return rho_m_z(z)/rho_c_z(z)

# Virial overdensity
def Deltac(z):
    d = Omz(z) -1.
    return 18.*np.pi**2. + 82.*d - 39.*d**2.

# Virial radius, kpc/h, M in Msun/h, physical units
def Rvir(M,z):
    return 0.784*(M/1.e8)**(1./3.)*(Omega_m/Omz(z)*Deltac(z)/(18.*np.pi**2.))**(-1./3.)*((1.+z)/10.)**(-1.)

# Circular velocity, km/s, Vcirc = sqrt(G*M/rvir)
def Vcirc(M,z):
    return 23.4*(M/1.e8)**(1./3.)*(Omega_m/Omz(z)*Deltac(z)/(18.*np.pi**2.))**(1./6.)*((1.+z)/10.)**(1./2)

# Virial temperature, K
# Barkana & Loeb 2001 has factor 1.98, while Shimabukuro et al 2014 has 1.32, as factor 3/2 lower
def Tvir(M,z):
    return 1.98e4*(mu/0.6)*(M/1.e8)**(2./3.)*(Omega_m*Deltac(z)/(Omz(z)*18.*np.pi**2.))**(1./3.)*((1.+z)/10.)

# Virial mass from a given virial temperature, Msun/h, Tvir in K
def Mmin(Tvir,z):
    return 1.e8*(Omega_m*Deltac(z)/(Omz(z)*18.*np.pi**2.))**(-1./2.)*(Tvir/1.98e4)**(3./2.)*((1.+z)/10.)**(-3./2.)

# Halo concentration parameter
def con(M,z,model="ishiyama21"):
    if model=="Bullock2001":      # Employs fit from Bullock et al. 2001, to match at z=0
        return 9./(1.+z)*(M/1.5e13)**(-0.13)
    else:     # Employs fit from colossus. Default here is taken as "ishiyama21"
        return concentration(M, "vir", z=z, model=model)

# Halo mass function in LCDM cosmology, dndlnM, units of (h/Mpc)^3, using the Sheth-Tormen parameterization
def dndlnM_CDM(M,z):
    return mass_function.massFunction(M,z, q_out = "dndlnM", model = "sheth99")

#--- DENSITY PROFILES ---#

# Integral appearing in NFW profile
def Fint(y):
    return np.log(1.+y)-y/(1.+y)

# Squared escape velocity, (km/s)^2
def Vesc2(M,z,y,r):
    x = r/Rvir(M,z)
    return 2.*Vcirc(M,z)**2.*( Fint(y*x) + y*x/(1.+y*x) )/(x*Fint(y))

# Central gas density, g/cm^3 . It fulfills M_gas = (Omega_b/Omega_m)*M_nfw
def rho_gas_central(z,y):
    A = 2.*y/Fint(y)
    lntvec = np.linspace(np.log(1.e-10),np.log(y))
    integral = integrate.simps((1.+np.exp(lntvec))**(A/np.exp(lntvec))*np.exp(lntvec)**3.,lntvec)
    return (Omega_b/Omega_m)*rho_c_z(z)*Deltac(z)/3.*y**3.*np.exp(A)/integral

# Gas density profile, obtained assuming an isothermal halo in hydrostatic equilibrium, g/cm^3
def rho_gas(M,z,y,r):
    Vesc2_0 = 2.*Vcirc(M,z)**2.*y/Fint(y)
    return rho_gas_central(z,y)*np.exp(-mu*m_p*(Vesc2_0-Vesc2(M,z,y,r))*1.e10/(2.*Tvir(M,z)*kelvintoerg))

# NFW density profile, g/cm^3
def rho_nfw(M,z,y,r):
    rho_0 = rho_c_z(z)*Deltac(z)/3.*y**3./Fint(y)
    x = y*r/Rvir(M,z)
    return rho_0/x/(1.+x)**2.



#--- 21 CM ---#

# Scattering rate for Hydrogen-Hydrogen collisions
def kappa_HH(T):
    if T>10:
        kh = 3.1e-11*T**(0.357)*np.exp(-32./T)
    elif T<=10 and T>1:
        kh = 3.6e-16*T**(3.640)*np.exp(6.035/T)
    else:
        kh = 3.6e-16*1.**(3.640)*np.exp(6.035/1.)
    return kh

# Collision coupling coefficient y
def yColl(z,nH,T):
    return Tstar/T/A10*kappa_HH(T)*nH

# Spin temperature, K
def Ts(z,nH,Tk):
    yc = yColl(z,nH,Tk)
    return (Tcmb0*(1.+z) + yc*Tk)/(1. + yc)

# 21 cm optical depth
def OptDepth21(M,z,y,al):
    prefactor = 3.*hplanck*c**3.*A10/(32.*np.pi*nu21**2.)   # cm^3*erg/s
    RRmax = np.sqrt(1.-(al/Rvir(M,z))**2.)        # Maximum radius to integrate
    if RRmax==0.:   # avoid division by 0
        return 0.
    logRR = np.linspace(np.log(1.e-3),np.log(RRmax),100)
    RR = np.exp(logRR)
    rvec = np.sqrt(RR**2. + (al/Rvir(M,z))**2.)
    Tk = Tvir(M,z)
    nH = rho_gas(M,z,y,rvec*Rvir(M,z))*(1.-Y_He)/m_p
    Tspin = Ts(z,nH,Tk)
    b = np.sqrt(2.*Tk*kelvintoerg/m_p)      # cm/s
    nu = nu21/(1.+z)
    expdoppler = 1./np.sqrt(np.pi)/b        # Doppler exponential factor
    #expdoppler = 1./np.sqrt(np.pi)/b*np.exp(-(c*(nu-nu21)/nu21)**2./b**2.)
    integral = 2.*Rvir(M,z)*integrate.simps( RR*nH/Tspin*expdoppler, logRR) # 2 is 'cause goes from 0 to Rmax
    #integral = 2.*Rvir(M,z)*integrate.quad(lambda logR: np.exp(logR)*rho_gas(M,z,y,np.sqrt(np.exp(logR)*2.+(al/Rvir(M,z))**2.)*Rvir(M,z))*(1.-Y_He)/m_p/Ts(z,rho_gas(M,z,y,np.sqrt(np.exp(logR)*2.+(al/Rvir(M,z))**2.)*Rvir(M,z))*(1.-Y_He)/m_p,Tk)*expdoppler, logRR[0], logRR[-1])[0]
    return prefactor*integral/kelvintoerg*MpcToCm/1.e3/hlittle

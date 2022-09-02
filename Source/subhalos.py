""""
Subhalos module
Includes the relevant routines for computing the subhalo contribution to the optical depth
Pablo Villanueva Domingo, Kenji Kadota
Started: October 2021
Last update: September 2022
"""


from Source.cosmo import *

# Set random seed for Monte Carlo integration
np.random.seed(0)

# Number of samples drawn within the domain for the Monte Carlo integraton. The bigger the more accurate but slower
Nsamples = 10000

# Subhalo mass function
# Masses Msub and Mhost in Msun     CHECK h factors!
# submf: Subhalo mass function parameterization
# Choose between "Ando19" (1903.11427) and "Dooley17" (1610.00708)
def subhalo_mf(Msub, z, Mhost, submf = "Ando19"):

    # Eq. 2 1610.00708
    if submf=="Dooley17":
        K0, chi = 1.88e-3, 1.87 # K0 in Msun^-1
        return K0*(Msub)**(-chi)*Mhost
    # Appendix A of 1903.11427
    elif submf=="Ando19":
        a_h = 6.0e-4*(1.+z)*Mhost*(2.-0.08*np.log10(Mhost))**2.*(np.log10(Mhost/1.e-5))**(-0.1)
        b_h = 8.0e-5*Mhost*(2.-0.08*np.log10(Mhost))*(np.log10(Mhost/1.e-5))**(-0.08*z)*(np.log10(Mhost/1.e-8))**(-1.)*(np.log10(Mhost/1.e18))**(2.)
        mc = 0.05*(1.+z)*Mhost
        alpha_h = 0.2 + 0.02*z
        beta_h = 3.
        return (a_h + b_h*Msub**alpha_h)/Msub**2.*np.exp(-(Msub/mc)**beta_h)

# Normalized density profile
def w_dens(Mhost,z,y,r,sub_distr):
    if sub_distr=="nfw":
        return rho_nfw(Mhost,z,y,r)/Mhost
    elif sub_distr=="uni":
        return rho_c_z(z)*Deltac(z)/Mhost

# Number of subhalos within a cylinder at a given impact parameter from the center of the host halo, using Monte Carlo integration
# Msub: host halo mass
# z: redshift
# y: concentration of the host halo
# al: impact parameter of the host halo (in kpc/h)
# Mhost: host halo mass
def num_subhalos_alpha(Msub,z,y,al,Mhost,sub_distr):

    Rhost, Rsub = Rvir(Mhost,z), Rvir(Msub,z)
    alplus, almin = min(al+Rsub,Rhost), max(al-Rsub,0.)
    unitfact = (1.e-3*MpcToCm)**3./Msun

    if Rhost < (al + Rsub):
        Nsub = 0.0

    else:

        zzmax = np.sqrt(Rhost ** 2 - (al + Rsub) ** 2)
        rrmax = Rsub
        thetamax = 2.0 * np.pi
        zzmin = 0.0
        rrmin = 0.0
        thetamin = 0.0

        domain = (thetamax - thetamin) * (zzmax - zzmin) * (rrmax - rrmin)

        rr = rrmax*np.sqrt(np.random.uniform(low=0., high=1., size=Nsamples))  # N values uniformly drawn from a to b
        theta = np.random.uniform(low=thetamin, high=thetamax, size=Nsamples)  # N values uniformly drawn from c to d
        zz = np.random.uniform(low=zzmin, high=zzmax, size=Nsamples)  # N values uniformly drawn from 0 to 2Pi

        rad = np.sqrt(rr ** 2 * np.sin(theta) ** 2 + (al + rr * np.cos(theta)) ** 2 + zz ** 2)

        integrand = w_dens(Mhost,z,y,rad,sub_distr)*rr

        Nsub = 2 * (np.sum(integrand) * domain) / Nsamples * subhalo_mf(Msub, z, Mhost) * unitfact

    return Nsub

num_subhalos_alpha = np.vectorize(num_subhalos_alpha)

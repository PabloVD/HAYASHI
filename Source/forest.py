""""
Main 21 cm forest class
Author: Pablo Villanueva Domingo
Last update: September 2022
"""

import numpy as np
import os
from tqdm import tqdm
from scipy import integrate, interpolate

from Source.constants import *
from Source import cosmo
from Source import subhalos
from Source import lss

#---------------------------
# 21cm Forest general class
#---------------------------
class Forest():
    def __init__(self,
                 z,                         # z: redshift of evaluation
                 Tk,                        # Tk: gas temperature of the IGM at that redshift
                 use_subhalos = True,       # use_subhalos: 1 for including the subhalo contribution, 0 otherwise
                 dndlnM = lss.dndlnM_CDM, # dndlnM: halo mass function for host halos (default: Sheth-Thormen)
                 sub_distr = "nfw",         # sub_distr: distribution of subhalos within the halo. Choose between "nfw" for the NFW profile, and "uni" for uniform distribution
                 tidal = False             # tidal: if True, use tidal disruption in outer parts of subhalos
                 ):

        # Redshift
        self.z = z
        # Temperature of the IGM at redshift z
        self.Tk = Tk
        # Flag, =True to take into account subhalo contribution
        self.use_subhalos = use_subhalos
        # Halo mass function
        self.dndlnM = dndlnM
        # Subhalo distribution within the host halo
        self.sub_distr = sub_distr
        # Consider tidal disruption in outer parts of subhalos
        self.tidal = tidal

        #--- Number of bins ---#
        # Number of bins of halo masses
        self.nummass = 60
        # Number of bins of the impact parameter
        self.numalf = 60
        # Number of bins of optical depth values
        self.numtau = 200

        # Jeans mass
        self.MJ = cosmo.MJeans(self.z,self.Tk)
        # Maximum mass for minihalos, taken as the minimum mass to host star formation, corresponding to a given virial temperature of 10^4 K
        self.Mmax = cosmo.Mmin(1.e4,z)

        # Mass, impact parameter and optical depth arrays
        self.logMvec = np.linspace(np.log(cosmo.MJeans(z,cosmo.Tk_ad(z))),np.log(self.Mmax),num=self.nummass)
        #self.logimpparam = np.linspace(np.log(1.e-5),0.,num=self.numalf)
        self.Mvec = np.exp(self.logMvec)
        #self.impparam = np.exp(self.logimpparam)   # impact parameter normalized to virial radius
        self.impparam = np.linspace(0.,1.,num=self.numalf)
        self.tauvec = np.logspace(-3,0,num=self.numtau)

        # Tidal radius: td_fac*r_s
        self.td_fac = 0.77  # value from astro-ph/0203004

        # Optical depth array (shape: tau[mass, impact param])
        self.tau_host_file = "Data/tau_z_{:.1f}".format(z)+"_nummass_"+str(self.nummass)+"_numalf_"+str(self.numalf)+".npy"
        if os.path.exists(self.tau_host_file):
            self.tau_host_arr = np.load(self.tau_host_file)
        else:
            print("Optical depth file for z="+str(z)+" does not exist. Generating...")
            self.tau_host_arr = self.get_optical_depth()

        # Subhalo contribution to optical depth array (shape: tau[mass, impact param])
        if self.use_subhalos:
            self.tau_sub_file = "Data/tau_sub_z_{:.1f}_Tk_{:.3e}".format(z,Tk)+"_"+self.sub_distr+"_tidal_"+str(int(self.tidal))+"_nummass_"+str(self.nummass)+"_numalf_"+str(self.numalf)+".npy"
            if os.path.exists(self.tau_sub_file):
                self.tau_sub_arr = np.load(self.tau_sub_file)
            else:
                print("Optical depth subhalo contribution file for z="+str(z)+", Tk={:.3e}".format(Tk)+" K does not exist. Generating...")
                self.tau_sub_arr = self.get_optical_depth_subhalos()
        else:
            self.tau_sub_arr = np.zeros((self.nummass,self.numalf))

        # Total optical depth
        self.tau_tot = self.tau_host_arr + self.tau_sub_arr

        # Boost factor for the optical depth
        self.boost_factor_arr = np.where(self.tau_host_arr==0., 0., self.tau_sub_arr/self.tau_host_arr)

        # Maximum impact parameter matrix (shape: tau[tau, mass])
        self.max_impact_arr = self.max_impact_param()

        # Create 2d interpolations of the optical depths, boost factor and maximum impact parameter
        # ..._int_logM expects the logarithm of mass as input, while ..._int expects as input the mass
        self.tau_host_int_logM = interpolate.interp2d(self.logMvec, self.impparam, np.transpose(self.tau_host_arr),kind="linear")
        self.tau_host_int = lambda M, al: self.tau_host_int_logM(np.log(M),al)
        self.tau_sub_int_logM = interpolate.interp2d(self.logMvec, self.impparam, np.transpose(self.tau_sub_arr),kind="linear")
        self.tau_sub_int = lambda M, al: self.tau_sub_int_logM(np.log(M),al)
        self.tau_tot_int = lambda M, al: self.tau_host_int(M, al) + self.tau_sub_int(M, al)
        self.boost_factor_int_logM = interpolate.interp2d(self.logMvec, self.impparam, np.transpose(self.boost_factor_arr),kind="linear")
        self.boost_factor_int = lambda M, al: self.boost_factor_int_logM(np.log(M),al)
        self.max_impact_int_log = interpolate.interp2d(np.log(self.tauvec), self.logMvec, np.transpose(self.max_impact_arr),kind="linear")
        self.max_impact_int = lambda tau, M: self.max_impact_int_log(np.log(tau),np.log(M))

    # Method to compute the optical depth array
    def get_optical_depth(self):

        tau_arr = np.zeros((self.nummass,self.numalf))

        pbar = tqdm(self.Mvec, total=self.Mvec.shape[0], position=0, leave=True, desc=f"Generating optical depth file")

        for iM, M in enumerate(pbar):
            for ia, impa in enumerate(self.impparam):

                tau_arr[iM,ia] = cosmo.OptDepth21(M,self.z,cosmo.con(M,self.z),impa*cosmo.Rvir(M,self.z))

        np.save(self.tau_host_file,tau_arr)

        return tau_arr

    # Method to compute the subhalo contribution to the optical depth array
    def get_optical_depth_subhalos(self):

        # Average tau
        tau_av = np.zeros(self.nummass)
        for iM, Msub in enumerate(self.Mvec):
            if self.tidal:
                stepfactor = np.heaviside(self.td_fac/cosmo.con(Msub,self.z) - self.impparam,1.)    # Max alpha: td_fac*rs (impparam is given as impact param normalized to Rvir, so a concentration factor appears)
            else:
                stepfactor = 1.
            tau_av[iM] = 2.*integrate.simps(self.tau_host_arr[iM,:]*self.impparam*stepfactor, self.impparam)

        tau_sub_arr = np.zeros((self.nummass,self.numalf))

        pbar = tqdm(self.Mvec, total=self.Mvec.shape[0], position=0, leave=True, desc=f"Generating optical depth subhalo contribution file")

        # For each host halo mass M
        for iM, Mhost in enumerate(pbar):
            if iM>0:
                for ia, impa in enumerate(self.impparam):

                    # Subhalo masses are those halo masses below the host halo mass
                    Msubs = self.Mvec[:iM]
                    avopdeps = tau_av[:iM]
                    # Number of subhalos array for each subhalo mass
                    numsubsalp = subhalos.num_subhalos_alpha(Msubs,self.z,cosmo.con(Mhost,self.z),impa*cosmo.Rvir(Mhost,self.z),Mhost,self.sub_distr)
                    #numsubsalp = np.array( [ subhalos.num_subhalos_alpha(Msub,self.z,cosmo.con(Mhost,self.z),impa*cosmo.Rvir(Mhost,self.z),Mhost,self.sub_distr) for Msub in Msubs ] )
                    # Write Jeans mass lower bound as a heaviside step function
                    stepfactor = np.heaviside(Msubs - self.MJ,1.)
                    tau_sub_arr[iM, ia] = integrate.simps(avopdeps*numsubsalp*stepfactor, Msubs)

        np.save(self.tau_sub_file,tau_sub_arr)

        return tau_sub_arr

    # Compute the maximum impact parameter for a given optical depth and mass
    def max_impact_param(self):

        max_impact_arr = np.zeros((self.numtau,self.nummass))

        for it, tau in enumerate(self.tauvec):
            maximpactvec = []
            for iM, M in enumerate(self.Mvec):
                maximpact = 0.
                for ia, impa in enumerate(self.impparam):
                    opdep = self.tau_tot[iM,ia]
                    if opdep>=tau:
                        maximpact = impa
                max_impact_arr[it, iM] = cosmo.Rvir(M,self.z)*maximpact

        return np.array(max_impact_arr)

    # Get the cumulative number of absorbers dN/dz and its logarithmic derivative taud^2N/dz/dtau
    def num_absorbers(self):

        cumulvec = np.zeros(self.numtau)
        prefactor = (1.+self.z)**2.*cosmo.drdz(self.z)/1.e6/(MpcToCm/hlittle)
        dlogtau = np.log(self.tauvec[1])-np.log(self.tauvec[0])

        for tt, tau in enumerate(self.tauvec):
            maximpactvec = self.max_impact_arr[tt]

            stepfactor = np.heaviside(self.logMvec - np.log(self.MJ),1.)
            intcumul = integrate.simps( self.dndlnM(np.exp(self.logMvec),self.z)*np.pi*maximpactvec**2.*stepfactor, self.logMvec )

            cumulvec[tt] = intcumul*prefactor

        # Compute the derivative of the cumulative number of absorbers with respect to tau
        dercumulvec = (cumulvec[1:] - cumulvec[:-1])/dlogtau

        self.cumulvec = cumulvec
        self.dercumulvec = dercumulvec

        return cumulvec, dercumulvec

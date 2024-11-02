""""
Utility functions
Pablo Villanueva Domingo
Started: October 2018
Last update: July 2022
"""

import numpy as np
import math
from scipy import interpolate
from hayashi.cosmo import Tk_ad
import os

# Create some folders for outputs if they are absent
for path in ["Plots", "Data"]:
        if not os.path.exists(path):
            os.mkdir(path)


#--- MISC ---#

# Write a number in Latex scientific notation
def scinot(x):
    exp = int(math.floor(math.log10(abs(x))))
    prefactor = x / 10**exp
    if exp==0:
        return r"{:.0f}".format(prefactor)
    if exp==1:
        return r"{:.0f}".format(prefactor*10.)
    elif prefactor == 1.:
        return r"$10^{"+str(exp)+"}$"
    else:
        return r"${:.0f}".format(prefactor)+" \\times 10^{"+str(exp)+"}$"


# Import the IGM temperature from interpolation
def GetTemperature(z, zetax):

    if zetax==0.:
        return Tk_ad(z)

    # Take temperature from file if it exists
    ThermalHisF = 'ThermalHistory_Xi_30_Tmin_1.000e+04_Rfmp_15_chiUV_{:.2e}_Nalpha_4.00e+03.dat'.format(zetax)
    file = np.loadtxt("21cmFiles/"+ThermalHisF,unpack=True)
    Temp = interpolate.interp1d(file[0], file[2], fill_value="extrapolate")
    Tk = Temp(z)
    return Tk

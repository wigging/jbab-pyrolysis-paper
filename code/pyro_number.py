"""
here
"""

import chemics as cm
import numpy as np

# Parameters
# ----------------------------------------------------------------------------

from params import dp_feed
from params import k_feed
from params import rhop_feed

from params import dp_bed
from params import ep
from params import phi_bed
from params import rhop_bed

from params import press
from params import temp


cp_feed = 103.1 + (3.867 * temp)


# Functions
# ----------------------------------------------------------------------------

def biot(h, r, k):
    """
    Calculate the dimensionless Biot number represented by Bi.

    Parameters
    ----------
    h : float
        Convective heat transfer coefficient [W/m²K]
    r : float
        Radius of the biomass particle [m]
    k : float
        Thermal conductivity of the biomass particle [W/mK]

    Returns
    -------
    Bi : float
        Biot number [-]
    """
    Bi = (h * r) / k
    return Bi


def pyro1(k, kr, rho, cp, r):
    """
    Calculate the pyrolysis number Py I.

    Parameters
    ----------
    k : float
        Thermal conductivity of the biomass particle [W/mK]
    kr : float
        Rate constant [1/s]
    rho : float
        Density of the biomass particle [kg/m³]
    cp : float
        Heat capacity of the biomass particle [J/kgK]
    r : float
        Radius of the biomass particle [m]

    Returns
    -------
    py : float
        Pyrolysis number Py I [-]
    """
    py = k / (kr * rho * cp * (r**2))
    return py


def pyro2(h, kr, rho, cp, r):
    """
    Calculate the pyrolysis number Py II.

    Parameters
    ----------
    h : float
        Convective heat transfer coefficient [W/m²K]
    kr : float
        Rate constant [1/s]
    rho : float
        Density of the biomass particle [kg/m³]
    cp : float
        Heat capacity of the biomass particle [J/kgK]
    r : float
        Radius of the biomass particle [m]

    Returns
    -------
    py : float
        Pyrolysis number Py II [-]
    """
    py = h / (kr * rho * cp * r)
    return py


# Diameter and radius of biomass particle
# ----------------------------------------------------------------------------

dps = []
wts = []

for x in dp_feed:
    dps.append(x['d'])
    wts.append(x['mf'])

# biomass particle average Sauter mean diameter [μm]
d_avg = np.average(dps, weights=wts)

# average biomass particle diameter in meters [m]
dp_avg = d_avg / 1e6

# radius of the biomass particle [m]
r = dp_avg / 2

# Reaction rate
# ----------------------------------------------------------------------------

# universal gas constant [kJ/(mol K)]
rconst = 0.008314

# biomass -> gas
a1 = 4.38e9
e1 = 152.7
k1 = a1 * np.exp(-e1 / (rconst * temp))

# biomass -> char
a2 = 3.27e6
e2 = 111.7
k2 = a2 * np.exp(-e2 / (rconst * temp))

# biomass -> tar
a3 = 1.08e10
e3 = 148.0
k3 = a3 * np.exp(-e3 / (rconst * temp))

# overall reaction rate constant for biomass conversion [1/s]
kr = k1 + k2 + k3

# Biot and pyrolysis numbers
# ----------------------------------------------------------------------------

# gases for calculations
gas = ['N2', 'H2', 'H2O', 'CO', 'CO2', 'CH4']

# Umf, Reynolds number, Nusselt number, and convective heat transfer
# coefficient (h) for each gas
umf = []
reynolds = []
nusselt = []
hconv = []

bi = []
py1 = []
py2 = []

for g in gas:
    mw = cm.mw(g)
    rho_gas = cm.rhog(mw, press, temp)
    mu_gas = cm.mu_gas(g, temp) / 1e7    # convert µP to kg/(ms)

    umf_ergun = cm.umf_ergun(dp_bed, ep, mu_gas, phi_bed, rho_gas, rhop_bed)
    umf_grace = cm.umf_coeff(dp_bed, mu_gas, rho_gas, rhop_bed, coeff='grace')
    umf_rich = cm.umf_coeff(dp_bed, mu_gas, rho_gas, rhop_bed, coeff='rich')
    umf_wenyu = cm.umf_coeff(dp_bed, mu_gas, rho_gas, rhop_bed, coeff='wenyu')
    umf_avg = (umf_ergun + umf_grace + umf_rich + umf_wenyu) / 4
    umf.append(umf_avg)

    re = (rho_gas * umf_avg * dp_avg) / mu_gas
    reynolds.append(re)

    nu = 2 + (0.9 * re**0.62) * ((dp_avg / dp_bed)**0.2)
    nusselt.append(nu)

    if g == 'CH4':
        k_gas = cm.k_gas_organic(g, temp)
    else:
        k_gas = cm.k_gas_inorganic(g, temp)
    h = (k_gas * nu) / dp_avg
    hconv.append(h)

    biotnum = biot(h, r, k_feed)
    bi.append(biotnum)

    py1num = pyro1(k_feed, kr, rhop_feed, cp_feed, r)
    py1.append(py1num)

    py2num = pyro2(h, kr, rhop_feed, cp_feed, r)
    py2.append(py2num)

# Print
# ----------------------------------------------------------------------------

print('')
print(f'{"gas":8} {"Umf":8} {"Re":8} {"Nu":8} {"h":8}')
for i in range(len(gas)):
    print(f'{gas[i]:<8} {umf[i]:<8.2f} {reynolds[i]:<8.2f} {nusselt[i]:<8.2f} {hconv[i]:<8.2f}')

print('')
print(f'{"gas":8} {"Bi":8} {"Py I":8} {"Py II":8}')
for i in range(len(gas)):
    print(f'{gas[i]:<8} {bi[i]:<8.2f} {py1[i]:<8.2f} {py2[i]:<8.2f}')

# Plot
# ----------------------------------------------------------------------------

"""
Compare the Biot and pyrolysis numbers for different fluidization gases.
"""

import chemics as cm
import matplotlib.pyplot as plt
import numpy as np
import params as pm
from funcs import biot, pyro1, pyro2

# Parameters
# ----------------------------------------------------------------------------

# heat capacity of the biomass [J/(kg⋅K)]
cp_feed = 103.1 + (3.867 * pm.temp)

# gases for calculations
gas = ['N2', 'H2', 'H2O', 'CO', 'CO2', 'CH4']

# range of biomass particle diameters [µm] for calculations
dmin = 10
dmax = 5000
diams = np.linspace(dmin, dmax, 10)

# Diameter of biomass particle (d_avg)
# ----------------------------------------------------------------------------

dps = []
wts = []

for x in pm.dp_feed:
    dps.append(x['d'])
    wts.append(x['mf'])

# average biomass particle diameter [m]
d_avg = np.average(dps, weights=wts) / 1e6

# Reaction rate constant (kr)
# ----------------------------------------------------------------------------

# universal gas constant [kJ/(mol K)]
rconst = 0.008314

# biomass -> gas
a1 = 4.38e9
e1 = 152.7
k1 = a1 * np.exp(-e1 / (rconst * pm.temp))

# biomass -> char
a2 = 3.27e6
e2 = 111.7
k2 = a2 * np.exp(-e2 / (rconst * pm.temp))

# biomass -> tar
a3 = 1.08e10
e3 = 148.0
k3 = a3 * np.exp(-e3 / (rconst * pm.temp))

# overall reaction rate constant for biomass conversion [1/s]
kr = k1 + k2 + k3

# Biot and pyrolysis numbers relevant to each gas (Bi, PyI, PyII)
# ----------------------------------------------------------------------------

# store properties
h_gas = []
bi_gas = []
py1_gas = []
py2_gas = []

for g in gas:
    mw = cm.mw(g)
    rho_gas = cm.rhog(mw, pm.press, pm.temp)
    mu_gas = cm.mu_gas(g, pm.temp) / 1e7    # convert µP to kg/(ms)

    umf_ergun = cm.umf_ergun(pm.dp_bed, pm.ep, mu_gas, pm.phi_bed, rho_gas, pm.rhop_bed)
    umf_grace = cm.umf_coeff(pm.dp_bed, mu_gas, rho_gas, pm.rhop_bed, coeff='grace')
    umf_rich = cm.umf_coeff(pm.dp_bed, mu_gas, rho_gas, pm.rhop_bed, coeff='rich')
    umf_wenyu = cm.umf_coeff(pm.dp_bed, mu_gas, rho_gas, pm.rhop_bed, coeff='wenyu')
    umf_avg = (umf_ergun + umf_grace + umf_rich + umf_wenyu) / 4

    re = (rho_gas * umf_avg * d_avg) / mu_gas
    nu = 2 + (0.9 * re**0.62) * ((d_avg / pm.dp_bed)**0.2)

    if g == 'CH4':
        k_gas = cm.k_gas_organic(g, pm.temp)
    else:
        k_gas = cm.k_gas_inorganic(g, pm.temp)
    h = (k_gas * nu) / d_avg
    h_gas.append(h)

    bi = biot(h, d_avg / 2, pm.k_feed)
    bi_gas.append(bi)

    py1 = pyro1(pm.k_feed, kr, pm.rhop_feed, cp_feed, d_avg / 2)
    py1_gas.append(py1)

    py2num = pyro2(h, kr, pm.rhop_feed, cp_feed, d_avg / 2)
    py2_gas.append(py2num)


# Biot and pyrolysis numbers for range of diameters (Bi, PyI, PyII)
# ----------------------------------------------------------------------------


class Feedstock:

    def __init__(self, cp, d, h, kr, params, diams=0):
        self.d = d
        self.h = h
        self.kr = kr
        self.cp = cp
        self.k = params.k_feed
        self.rho = params.rhop_feed
        self.diams = diams

    @property
    def biot(self):
        r = self.d / 2
        bi = (self.h * r) / self.k
        return bi

    @property
    def pyro(self):
        r = self.d / 2
        pyI = self.k / (self.kr * self.rho * self.cp * (r**2))
        pyII = self.h / (self.kr * self.rho * self.cp * r)
        if self.biot < 1.0:
            return pyII
        else:
            return pyI

    def calc_from_diams(self):
        biot_diams = []
        pyro_diams = []

        for d in self.diams:
            d = d / 1e6
            r = d / 2
            bi = (self.h * r) / self.k
            biot_diams.append(bi)

            pyI = self.k / (self.kr * self.rho * self.cp * (r**2))
            pyII = self.h / (self.kr * self.rho * self.cp * r)
            if bi < 1.0:
                pyro_diams.append(pyII)
            else:
                pyro_diams.append(pyI)

        return biot_diams, pyro_diams


feed_n2 = Feedstock(cp=cp_feed, d=d_avg, h=370, kr=kr, params=pm, diams=diams)
feed_h2 = Feedstock(cp=cp_feed, d=d_avg, h=2200, kr=kr, params=pm, diams=diams)

biot_n2, pyro_n2 = feed_n2.calc_from_diams()
biot_h2, pyro_h2 = feed_h2.calc_from_diams()

# Print
# ----------------------------------------------------------------------------

print(f'{"gas":8} {"Bi":8} {"Py I":8} {"Py II":8}')
for i in range(len(gas)):
    print(f'{gas[i]:<8} {bi_gas[i]:<8.2f} {py1_gas[i]:<8.2f} {py2_gas[i]:<8.2f}')

# Plot
# ----------------------------------------------------------------------------

# Figure 1
fig, ax = plt.subplots(tight_layout=True)
for i in range(len(gas)):
    bi = bi_gas[i]
    pyI = py1_gas[i]
    pyII = py2_gas[i]
    if bi < 1.0:
        ax.plot(bi, pyII, 'o', label=gas[i])
    else:
        ax.plot(bi, pyI, '^', label=gas[i])
ax.set_xlabel('Biot Number, Bi [-]')
ax.set_ylabel('Pyrolysis Number, Py [-]')

ax.text(0.2, 0.91, 'kinetics limited\nisothermal', ha='center', transform=ax.transAxes)
ax.text(0.8, 0.91, 'kinetics limited\nnon-isothermal', ha='center', transform=ax.transAxes)
ax.text(0.2, 0.03, 'convection limited', ha='center', transform=ax.transAxes)
ax.text(0.8, 0.03, 'conduction limited', ha='center', transform=ax.transAxes)
ax.axvline(1, c='k', ls='-.')
ax.axvspan(10**-1, 10**1, color='0.9')
ax.axhline(1, c='k', ls='-.')
ax.axhspan(10**-1, 10**1, color='0.9')
ax.grid(color='0.9')
ax.legend(loc='right', frameon=False)
ax.set_frame_on(False)
ax.set_xlim(10**-4, 10**4)
ax.set_ylim(10**-4, 10**4)
ax.set_xscale('log')
ax.set_yscale('log')
ax.tick_params(color='0.9')
plt.minorticks_off()

# Figure 2
fig, ax = plt.subplots(tight_layout=True)
ax.plot(biot_n2, pyro_n2, marker='.')
ax.plot(biot_h2, pyro_h2, marker='.')
ax.plot(feed_n2.biot, feed_n2.pyro, 'k^')
ax.plot(feed_h2.biot, feed_h2.pyro, 'k^')
ax.set_xlabel('Biot number, Bi [-]')
ax.set_ylabel('Pyrolysis number, Py [-]')

ax.text(biot_n2[0] - 0.008, pyro_n2[0], 'N₂')
ax.text(biot_h2[0] - 0.05, pyro_h2[0], 'H₂')

ax.text(biot_h2[0] + 0.02, pyro_h2[0], f'{diams[0]:.2g} µm')
ax.text(biot_h2[-1] + 10, pyro_h2[-1], f'{diams[-1] / 1000:.1g} mm')

ax.text(0.2, 0.91, 'kinetics limited\nisothermal', ha='center', transform=ax.transAxes)
ax.text(0.8, 0.91, 'kinetics limited\nnon-isothermal', ha='center', transform=ax.transAxes)
ax.text(0.2, 0.03, 'convection limited', ha='center', transform=ax.transAxes)
ax.text(0.8, 0.03, 'conduction limited', ha='center', transform=ax.transAxes)
ax.axvline(1, c='k', ls='-.')
ax.axvspan(10**-1, 10**1, color='0.9')
ax.axhline(1, c='k', ls='-.')
ax.axhspan(10**-1, 10**1, color='0.9')
ax.grid(color='0.9')
ax.set_frame_on(False)
ax.set_xlim(10**-4, 10**4)
ax.set_ylim(10**-4, 10**4)
ax.set_xscale('log')
ax.set_yscale('log')
ax.tick_params(color='0.9')
plt.minorticks_off()

plt.show()

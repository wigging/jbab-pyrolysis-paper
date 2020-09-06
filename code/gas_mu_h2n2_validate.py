"""
Compare H₂/N₂ gas viscosity data to calculated values. Note that P is units of
poise.
"""

import chemics as cm
import matplotlib.pyplot as plt
import numpy as np

from funcs.mu_brokaw import mu_brokaw
from funcs.mu_davidson import mu_davidson
from funcs.mu_wilke import mu_wilke

# Parameters
# ----------------------------------------------------------------------------

# H₂ and N₂ gas viscosity at 291.1 K (18°C) in units of P x 10^7
# values taken from data points (see below)
mu_h2 = 877
mu_n2 = 1752

# molecular weight [g/mol]
mw_h2 = 2.016
mw_n2 = 28.014

# Data
# ----------------------------------------------------------------------------

# H₂/N₂ viscosity data from Table 1 at 291.1 K (18°C) in Itterbeek 1947 paper
# mu is gas viscosity in P x 10^7
# x is H₂ mole fraction of the H₂/N₂ mixture
# y is H₂ mass fraction of the H₂/N₂ mixture
data_itterbeek = {
    'mu': [877, 1251, 1560, 1660, 1677, 1742, 1752],
    'x': [1, 0.84, 0.559, 0.38, 0.241, 0.134, 0],
    'y': []
}

# calculate H₂ mass fraction
for xh2 in data_itterbeek['x']:
    xn2 = 1 - xh2
    yh2 = cm.molefrac_to_massfrac([xh2, xn2], [mw_h2, mw_n2])[0]
    data_itterbeek['y'].append(yh2)

# H₂/N₂ viscosity data from Table 4 at 19°C in Trautz 1929 paper
data_trautz = {
    'mu': [874, 1305, 1472, 1598, 1703, 1739],
    'x': [1, 0.8077, 0.6672, 0.5053, 0.2021, 0],
    'y': []
}

# calculate H₂ mass fraction
for xh2 in data_trautz['x']:
    xn2 = 1 - xh2
    yh2 = cm.molefrac_to_massfrac([xh2, xn2], [mw_h2, mw_n2])[0]
    data_trautz['y'].append(yh2)

# Calculate gas mixture viscosity
# ----------------------------------------------------------------------------

# store calculated viscosity of gas mixture and associated mass fraction
mu_h2n2 = {
    'brokaw': [],
    'davidson': [],
    'graham': [],
    'herning': [],
    'wilke': []
}

y_h2 = []

# H₂ mole fractions for calculations
# endpoints chosen to avoid division by zero
x_h2 = np.linspace(0.0001, 0.9999)

for xh2 in x_h2:
    xn2 = 1.0 - xh2

    mu0 = mu_brokaw([mu_h2, mu_n2], [mw_h2, mw_n2], [xh2, xn2])
    mu_h2n2['brokaw'].append(mu0)

    mu1 = mu_davidson([mu_h2, mu_n2], [mw_h2, mw_n2], [xh2, xn2])
    mu_h2n2['davidson'].append(mu1)

    mu2 = cm.mu_graham([mu_h2, mu_n2], [xh2, xn2])
    mu_h2n2['graham'].append(mu2)

    mu3 = cm.mu_herning([mu_h2, mu_n2], [mw_h2, mw_n2], [xh2, xn2])
    mu_h2n2['herning'].append(mu3)

    mu4 = mu_wilke([mu_h2, mu_n2], [mw_h2, mw_n2], [xh2, xn2])
    mu_h2n2['wilke'].append(mu4)

    yh2 = cm.molefrac_to_massfrac([xh2, xn2], [mw_h2, mw_n2])[0]
    y_h2.append(yh2)

# Plot
# ----------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(figsize=(10, 4.8), nrows=1, ncols=2, sharey=True, tight_layout=True)

ax1.plot(data_itterbeek['x'], data_itterbeek['mu'], 'ko', label='Itterbeek')
ax1.plot(data_trautz['x'], data_trautz['mu'], 'k^', label='Trautz')
ax1.plot(x_h2, mu_h2n2['brokaw'], label='Brokaw')
ax1.plot(x_h2, mu_h2n2['davidson'], label='Davidson')
ax1.plot(x_h2, mu_h2n2['graham'], label='Graham')
ax1.plot(x_h2, mu_h2n2['herning'], label='Herning')
ax1.plot(x_h2, mu_h2n2['wilke'], label='Wilke')
ax1.set_xlabel('H₂ mole fraction [-]')
ax1.set_ylabel('Dynamic viscosity [P x 10$^7$]')
ax1.grid(color='0.9')
ax1.set_frame_on(False)
ax1.tick_params(color='0.9')

ax2.plot(data_itterbeek['y'], data_itterbeek['mu'], 'ko', label='Itterbeek')
ax2.plot(data_trautz['y'], data_trautz['mu'], 'k^', label='Trautz')
ax2.plot(y_h2, mu_h2n2['brokaw'], label='Brokaw')
ax2.plot(y_h2, mu_h2n2['davidson'], label='Davidson')
ax2.plot(y_h2, mu_h2n2['graham'], label='Graham')
ax2.plot(y_h2, mu_h2n2['herning'], label='Herning')
ax2.plot(y_h2, mu_h2n2['wilke'], label='Wilke')
ax2.set_xlabel('H₂ mass fraction [-]')
ax2.grid(color='0.9')
ax2.legend(loc='best', frameon=False)
ax2.set_frame_on(False)
ax2.tick_params(color='0.9')

fig.savefig('../tex/figures/gas-mu-h2n2-validate.pdf')

plt.show()

"""
Compare H₂/O₂ gas viscosity data to calculated values. Note that P is units of
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

# H₂ and O₂ gas viscosity at 293.6 K in units of P x 10^7
# values taken from data points (see below)
mu_h2 = 885
mu_o2 = 2040

# molecular weight [g/mol]
mw_h2 = 2.016
mw_o2 = 31.998

# Data
# ----------------------------------------------------------------------------

# H₂/O₂ mixture data from Table 2 at 293.6 K (20°C) in Itterbeek 1947 paper
# mu is gas viscosity in units of P x 10^7
# x is H₂ mole fraction of the H₂/O₂ mixture
# y is H₂ mass fraction of the H₂/O₂ mixture
data_itterbeek = {
    'mu': [885, 1409, 1615, 1739, 1868, 1954, 2040],
    'x': [1, 0.839, 0.727, 0.62, 0.473, 0.33, 0],
    'y': []
}

# calculate H₂ mass fraction
for xh2 in data_itterbeek['x']:
    xo2 = 1 - xh2
    yh2 = cm.molefrac_to_massfrac([xh2, xo2], [mw_h2, mw_o2])[0]
    data_itterbeek['y'].append(yh2)

# Calculate gas mixture viscosity
# ----------------------------------------------------------------------------

# store calculated viscosity of gas mixture and associated mass fraction
mu_h2o2 = {
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
    xo2 = 1 - xh2

    mu0 = mu_brokaw([mu_h2, mu_o2], [mw_h2, mw_o2], [xh2, xo2])
    mu_h2o2['brokaw'].append(mu0)

    mu1 = mu_davidson([mu_h2, mu_o2], [mw_h2, mw_o2], [xh2, xo2])
    mu_h2o2['davidson'].append(mu1)

    mu2 = cm.mu_graham([mu_h2, mu_o2], [xh2, xo2])
    mu_h2o2['graham'].append(mu2)

    mu3 = cm.mu_herning([mu_h2, mu_o2], [mw_h2, mw_o2], [xh2, xo2])
    mu_h2o2['herning'].append(mu3)

    mu4 = mu_wilke([mu_h2, mu_o2], [mw_h2, mw_o2], [xh2, xo2])
    mu_h2o2['wilke'].append(mu4)

    yh2 = cm.molefrac_to_massfrac([xh2, xo2], [mw_h2, mw_o2])[0]
    y_h2.append(yh2)

# Plot
# ----------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(figsize=(10, 4.8), nrows=1, ncols=2, sharey=True, tight_layout=True)

ax1.plot(data_itterbeek['x'], data_itterbeek['mu'], 'ko', label='Itterbeek')
ax1.plot(x_h2, mu_h2o2['brokaw'], label='Brokaw')
ax1.plot(x_h2, mu_h2o2['davidson'], label='Davidson')
ax1.plot(x_h2, mu_h2o2['graham'], label='Graham')
ax1.plot(x_h2, mu_h2o2['herning'], label='Herning')
ax1.plot(x_h2, mu_h2o2['wilke'], label='Wilke')
ax1.set_xlabel('H₂ mole fraction [-]')
ax1.set_ylabel('Dynamic viscosity [P x 10$^7$]')
ax1.grid(color='0.9')
ax1.set_frame_on(False)
ax1.tick_params(color='0.9')

ax2.plot(data_itterbeek['y'], data_itterbeek['mu'], 'ko', label='Itterbeek')
ax2.plot(y_h2, mu_h2o2['brokaw'], label='Brokaw')
ax2.plot(y_h2, mu_h2o2['davidson'], label='Davidson')
ax2.plot(y_h2, mu_h2o2['graham'], label='Graham')
ax2.plot(y_h2, mu_h2o2['herning'], label='Herning')
ax2.plot(y_h2, mu_h2o2['wilke'], label='Wilke')
ax2.set_xlabel('H₂ mass fraction [-]')
ax2.grid(color='0.9')
ax2.legend(loc='best', frameon=False)
ax2.set_frame_on(False)
ax2.tick_params(color='0.9')

fig.savefig('../tex/figures/gas-mu-h2o2-validate.pdf')

plt.show()

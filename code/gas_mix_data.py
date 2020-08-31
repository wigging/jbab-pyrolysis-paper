"""
Compare experiment data of gas mixture viscosity to calculated values. Note
that P is units of poise.
"""

import chemics as cm
import matplotlib.pyplot as plt
import numpy as np

from funcs import mu_davidson

# Parameters
# ----------------------------------------------------------------------------

# H₂ and N₂ gas viscosity at 291.1 K in units of P x 10^7
# O₂ gas viscosity at 293.6 K
# values from https://www.lmnoeng.com/Flow/GasViscosity.php
mu_h2 = 8.7027205e-5 * 1e7
mu_n2 = 1.7375524e-4 * 1e7
mu_o2 = 2.0254787e-4 * 1e7

# molecular weight [g/mol]
mw_h2 = 2.016
mw_n2 = 28.014
mw_o2 = 31.998

# Data
# ----------------------------------------------------------------------------

# H₂/N₂ mixture data from Table 1 at 291.1 K (18°C) in Itterbeek 1947 paper
# H₂/O₂ mixture data from Table 2 at 293.6 K (20°C) in Itterbeek 1947 paper
# mu is gas viscosity in P x 10^7
# h2 is percent of H₂ in mixture
data_itterbeek = {
    'h2n2': {
        'mu': np.array([877, 1251, 1560, 1660, 1677, 1742, 1752]),
        'h2': np.array([100, 84, 55.9, 38, 24.1, 13.4, 0])
    },
    'h2o2': {
        'mu': np.array([885, 1409, 1615, 1739, 1868, 1954, 2040]),
        'h2': np.array([100, 83.9, 72.7, 62, 47.3, 33, 0])
    }
}

# H₂/N₂ mixture data from Table 4 at 19°C in Trautz 1929 paper
data_trautz = {
    'h2n2': {
        'mu': np.array([874, 1305, 1472, 1598, 1703, 1739]),
        'h2': np.array([100, 80.77, 66.72, 50.53, 20.21, 0])
    }
}

# Calculate gas mixture viscosity
# ----------------------------------------------------------------------------

# fractions of gas mixture components
x_h2 = np.linspace(0, 1)
x_n2 = 1 - x_h2
x_o2 = 1 - x_h2

# store calculated gas viscosity of mixture
mu_h2n2_calc = {
    'davidson': [],
    'graham': [],
    'herning': []
}

mu_h2o2_calc = {
    'davidson': [],
    'graham': [],
    'herning': []
}

for i in range(len(x_h2)):
    xh2 = x_h2[i]
    xn2 = x_n2[i]
    xo2 = x_o2[i]

    mu1 = mu_davidson([mu_h2, mu_n2], [mw_h2, mw_n2], [xh2, xn2])
    mu_h2n2_calc['davidson'].append(mu1)

    mu2 = cm.mu_graham([mu_h2, mu_n2], [xh2, xn2])
    mu_h2n2_calc['graham'].append(mu2)

    mu3 = cm.mu_herning([mu_h2, mu_n2], [mw_h2, mw_n2], [xh2, xn2])
    mu_h2n2_calc['herning'].append(mu3)

    mu4 = mu_davidson([mu_h2, mu_o2], [mw_h2, mw_o2], [xh2, xo2])
    mu_h2o2_calc['davidson'].append(mu4)

    mu5 = cm.mu_graham([mu_h2, mu_o2], [xh2, xo2])
    mu_h2o2_calc['graham'].append(mu5)

    mu6 = cm.mu_herning([mu_h2, mu_o2], [mw_h2, mw_o2], [xh2, xo2])
    mu_h2o2_calc['herning'].append(mu6)

# Plot
# ----------------------------------------------------------------------------

fig, ax = plt.subplots(tight_layout=True)
ax.plot(data_itterbeek['h2n2']['h2'] / 100, data_itterbeek['h2n2']['mu'], 'ko', label='Itterbeek')
ax.plot(data_trautz['h2n2']['h2'] / 100, data_trautz['h2n2']['mu'], 'k^', label='Trautz')
ax.plot(x_h2, mu_h2n2_calc['davidson'], label='Davidson')
ax.plot(x_h2, mu_h2n2_calc['graham'], label='Graham')
ax.plot(x_h2, mu_h2n2_calc['herning'], label='Herning')
ax.set_xlabel('H₂ fraction in H₂/N₂ mixture [-]')
ax.set_ylabel('Viscosity [P x 10$^7$]')
ax.grid(color='0.9')
ax.legend(loc='best')
ax.set_frame_on(False)
ax.tick_params(color='0.9')

fig, ax = plt.subplots(tight_layout=True)
ax.plot(data_itterbeek['h2o2']['h2'] / 100, data_itterbeek['h2o2']['mu'], 'ko', label='Itterbeek')
ax.plot(x_h2, mu_h2o2_calc['davidson'], label='Davidson')
ax.plot(x_h2, mu_h2o2_calc['graham'], label='Graham')
ax.plot(x_h2, mu_h2o2_calc['herning'], label='Herning')
ax.set_xlabel('H₂ fraction in H₂/O₂ mixture [-]')
ax.set_ylabel('Viscosity [P x 10$^7$]')
ax.grid(color='0.9')
ax.legend(loc='best')
ax.set_frame_on(False)
ax.tick_params(color='0.9')

plt.show()

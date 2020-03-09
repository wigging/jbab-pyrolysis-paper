"""
Calculations for H2 and N2 gas and mixtures of the two gases.
"""

import chemics as cm
import matplotlib.pyplot as plt
import numpy as np

# Parameters
# ----------------------------------------------------------------------------

# Gas formulas for calculations
h2 = 'H2'
n2 = 'N2'

# Gas pressure [Pa] and temperature [K]
pgas = 101_325
tgas = 773.15

# Temperature range for calculations [K]
temp = np.linspace(723.15, 823.15, 20)

# Molecular Weight
# ----------------------------------------------------------------------------

mw_h2 = cm.mw(h2)
mw_n2 = cm.mw(n2)

# Gas Viscosity of H2 and N2 Mixture at Temperature (tgas)
# ----------------------------------------------------------------------------

mu_mix = np.array([cm.mu_gas('H2', tgas), cm.mu_gas('N2', tgas)])
mw_mix = np.array([mw_h2, mw_n2])
x_mix = np.array([0.85, 0.15])

# Gas Viscosity for Range of Temperatures
# ----------------------------------------------------------------------------

mu_h2 = []
mu_n2 = []
mu_h2n2_h = []
mu_h2n2_g = []

for t in temp:
    mu_h2.append(cm.mu_gas(h2, t))
    mu_n2.append(cm.mu_gas(n2, t))
    mu_i = np.array([cm.mu_gas('H2', t), cm.mu_gas('N2', t)])
    mu_h2n2_h.append(cm.mu_herning(mu_i, mw_mix, x_mix))
    mu_h2n2_g.append(cm.mu_graham(mu_i, x_mix))

# Print
# ----------------------------------------------------------------------------

print(f"""
Parameters
----------
P gas \t {pgas:,} Pa \t\t gas pressure
T gas \t {tgas} K ({tgas - 273.15}°C) \t gas temperature
T min \t {temp[0]} K ({temp[0] - 273.15}°C) \t min gas temperature
T max \t {temp[-1]} K ({temp[-1] - 273.15}°C) \t max gas temperature

Viscosity at T gas
------------------
mu_graham \t {cm.mu_graham(mu_mix, x_mix):.2f} µP
mu_herning \t {cm.mu_herning(mu_mix, mw_mix, x_mix):.2f} µP
""")

# Plot
# ----------------------------------------------------------------------------

fig, ax = plt.subplots(tight_layout=True)
ax.plot(temp, mu_h2, marker='.', label='H₂')
ax.plot(temp, mu_n2, marker='.', label='N₂')
ax.plot(temp, mu_h2n2_h, marker='.', label='0.85H₂ 0.15N₂ (h)')
ax.plot(temp, mu_h2n2_g, marker='.', label='0.85H₂ 0.15N₂ (g)')
ax.set_xlabel('Temperature [K]')
ax.set_ylabel('Viscosity [µP]')
ax.grid(color='0.9')
ax.legend(frameon=False, loc='best')
ax.set_frame_on(False)
ax.tick_params(color='0.9')

plt.show()

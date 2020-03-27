"""
Calculations for H2 and N2 gas and mixtures of the two gases.
"""

import chemics as cm
import matplotlib.pyplot as plt
import numpy as np

# Parameters
# ----------------------------------------------------------------------------

from params import temp
from params import temp_min
from params import temp_max

# Gas viscosity of H2 and N2 mixture at temperature
# ----------------------------------------------------------------------------

# molecular weight [g/mol]
mw_h2 = cm.mw('H2')
mw_n2 = cm.mw('N2')

mu_mix = np.array([cm.mu_gas('H2', temp), cm.mu_gas('N2', temp)])
mw_mix = np.array([mw_h2, mw_n2])
x_mix = np.array([0.85, 0.15])

# Gas viscosity of H2 and N2 mixture for range of temperatures
# ----------------------------------------------------------------------------

# temperature range for calculations [K]
temps = np.linspace(temp_min, temp_max, 20)

mu_h2 = []
mu_n2 = []
mu_h2n2_h = []
mu_h2n2_g = []

for tk in temps:
    mu_h2.append(cm.mu_gas('H2', tk))
    mu_n2.append(cm.mu_gas('N2', tk))
    mu_i = np.array([cm.mu_gas('H2', tk), cm.mu_gas('N2', tk)])
    mu_h2n2_h.append(cm.mu_herning(mu_i, mw_mix, x_mix))
    mu_h2n2_g.append(cm.mu_graham(mu_i, x_mix))

# Print
# ----------------------------------------------------------------------------

print(f"""
Parameters
----------
temp        {temp} K
temp_min    {temp_min} K
temp_max    {temp_max} K

Calculated viscosity at {temp} K
--------------------------------
mu_graham   {cm.mu_graham(mu_mix, x_mix):.2f} µP
mu_herning  {cm.mu_herning(mu_mix, mw_mix, x_mix):.2f} µP
""")

# Plot
# ----------------------------------------------------------------------------

fig, ax = plt.subplots(tight_layout=True)
ax.plot(temps, mu_h2, marker='.', label='H₂')
ax.plot(temps, mu_n2, marker='.', label='N₂')
ax.plot(temps, mu_h2n2_h, marker='.', label='0.85H₂ 0.15N₂ (Herning)')
ax.plot(temps, mu_h2n2_g, marker='.', label='0.85H₂ 0.15N₂ (Graham)')
ax.set_xlabel('Temperature [K]')
ax.set_ylabel('Viscosity [µP]')
ax.grid(color='0.9')
ax.legend(frameon=False, loc='best')
ax.set_frame_on(False)
ax.tick_params(color='0.9')

plt.show()

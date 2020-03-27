"""
Plot density of nitrogen and hydrogen gas for a range of pressures and
temperatures. Plot viscosity of nitrogen and hydrogen gas for a range of
temperatures.
"""

import chemics as cm
import matplotlib.pyplot as plt
import numpy as np

# Parameters
# ----------------------------------------------------------------------------

from params import temp
from params import press

# Gas Density at constant temperature
# ----------------------------------------------------------------------------

mw_h2 = cm.mw('H2')
mw_n2 = cm.mw('N2')

rho_p_h2 = []
rho_p_n2 = []

pressures = np.linspace(80_000, 160_000, 20)

for p in pressures:
    rho_p_h2.append(cm.rhog(mw_h2, p, temp))
    rho_p_n2.append(cm.rhog(mw_n2, p, temp))

# Gas density at constant pressure
# ----------------------------------------------------------------------------

rho_t_h2 = []
rho_t_n2 = []

temps = np.linspace(723, 823, 20)

for t in temps:
    rho_t_h2.append(cm.rhog(mw_h2, press, t))
    rho_t_n2.append(cm.rhog(mw_n2, press, t))

# Mesh grid to compare pressure, temperature, density
# ----------------------------------------------------------------------------

x, y = np.meshgrid(pressures, temps)
z = cm.rhog(mw_n2, x, y)

# Gas viscosity
# ----------------------------------------------------------------------------

mu_h2 = []
mu_n2 = []

for t in temps:
    mu_h2.append(cm.mu_gas('H2', t))
    mu_n2.append(cm.mu_gas('N2', t))

# Print
# ----------------------------------------------------------------------------

print(f"""
Parameters
----------
temp    {temp} K
press   {press:,} Pa
""")

# Plot
# ----------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(1, 2, tight_layout=True)
ax1.plot(pressures / 1000, rho_p_h2, marker='.', label='H₂')
ax1.plot(pressures / 1000, rho_p_n2, marker='.', label='N₂')
ax1.text(0.05, 0.82, f'T = {temp} K', transform=ax1.transAxes)
ax1.set_xlabel('Pressure [kPa]')
ax1.set_ylabel('Gas Density [kg/m³]')
ax1.grid(color='0.9')
ax1.legend(loc='best')
ax1.set_frame_on(False)
ax1.tick_params(color='0.9')
ax2.plot(temps, rho_t_h2, marker='.', label='H₂')
ax2.plot(temps, rho_t_n2, marker='.', label='N₂')
ax2.text(0.5, 0.39, f'P = {press / 1000} kPa', transform=ax2.transAxes)
ax2.set_xlabel('Temperature [K]')
ax2.set_ylabel('Gas Density [kg/m³]')
ax2.grid(color='0.9')
ax2.legend(loc='center right')
ax2.set_frame_on(False)
ax2.tick_params(color='0.9')

fig, ax = plt.subplots(tight_layout=True)
cs = ax.contourf(x / 1000, y, z)
ax.set_frame_on(False)
ax.set_title('Nitrogen Gas')
ax.set_xlabel('Pressure [kPa]')
ax.set_ylabel('Temperature [K]')
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Gas Density [kg/m³]')

fig, ax = plt.subplots(tight_layout=True)
ax.plot(temps, mu_h2, marker='.', label='H₂')
ax.plot(temps, mu_n2, marker='.', label='N₂')
ax.set_xlabel('Temperature [K]')
ax.set_ylabel('Viscosity [µP]')
ax.grid(color='0.9')
ax.legend(loc='best')
ax.set_frame_on(False)
ax.tick_params(color='0.9')

plt.show()

"""

Compare properties such as molecular weight, viscosity, thermal conductivity,
and density of individual gases.

"""

import chemics as cm
import matplotlib.pyplot as plt

# Parameters
# ----------------------------------------------------------------------------

# properties are calculated for each gas item
gas = ['CO', 'CO2', 'H2', 'H2O', 'N2']

# pressure [Pa] and temperature [K] of gas in the reactor
p_gas = 101_325
tk_gas = 773.15

# Gas Properties
# ----------------------------------------------------------------------------

mw = []     # store molecular weight of each gas
k = []      # store thermal conductivity of each gas
mu = []     # store viscosity of each gas
rho = []    # store density of each gas

for i in range(len(gas)):
    mw.append(cm.mw(gas[i]))
    k.append(cm.k_gas_inorganic(gas[i], tk_gas))
    mu.append(cm.mu_gas(gas[i], tk_gas))
    rho.append(cm.rhog(mw[i], p_gas, tk_gas))

# Print
# ----------------------------------------------------------------------------

print(f"""
Parameters
----------
P gas \t {p_gas:,} Pa \t\t gas pressure
T gas \t {tk_gas} K ({tk_gas - 273.15}°C) \t gas temperature
""")

# Plot
# ----------------------------------------------------------------------------

sub = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')

xticks = range(len(gas))
xlabels = [g.translate(sub) for g in gas]

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 6), tight_layout=True)

ax1.bar(xticks, mw, align='center', color='C0')
ax1.set_frame_on(False)
ax1.set_xticks(xticks)
ax1.set_xticklabels(xlabels)
ax1.set_ylabel('Molecular weight [g/mol]')

ax2.bar(xticks, mu, align='center', color='C4')
ax2.set_frame_on(False)
ax2.set_xticks(xticks)
ax2.set_xticklabels(xlabels)
ax2.set_ylabel('Viscosity [µP]')

ax3.bar(xticks, rho, align='center', color='C9')
ax3.set_frame_on(False)
ax3.set_xticks(xticks)
ax3.set_xticklabels(xlabels)
ax3.set_ylabel('Density [kg/m³]')

ax4.bar(xticks, k, align='center', color='C1')
ax4.set_frame_on(False)
ax4.set_xticks(xticks)
ax4.set_xticklabels(xlabels)
ax4.set_ylabel('Thermal conductivity [W/(m K)]')

plt.show()

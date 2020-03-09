"""
Compare properties such as molecular weight, viscosity, and density of
individual gases.
"""

import chemics as cm
import matplotlib.pyplot as plt

# Parameters
# ----------------------------------------------------------------------------

# Gas formulas where each item is a single gas.
# CH4 is methane, C2H6 is ethane, and C3H8 is propane.
gas = ['CO', 'CO2', 'H2', 'N2']

# Pressure of gas in the reactor [Pa]
pgas = 101_325

# Temperature of gas in the reactor [K]
tgas = 773.15

# Gas Properties
# ----------------------------------------------------------------------------

ngas = len(gas)
mw = []
mu = []
rho = []

for i in range(ngas):
    mw.append(cm.mw(gas[i]))
    mu.append(cm.mu_gas(gas[i], tgas))
    rho.append(cm.rhog(mw[i], pgas, tgas))

# Print
# ----------------------------------------------------------------------------

print(f"""
Parameters
----------
P gas \t {pgas:,} Pa \t\t gas pressure
T gas \t {tgas} K ({tgas - 273.15}°C) \t gas temperature
""")

# Plot
# ----------------------------------------------------------------------------

sub = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')

xticks = range(len(gas))
xlabels = [g.translate(sub) for g in gas]

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), tight_layout=True)

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

plt.show()

"""
Compare properties such as molecular weight, viscosity, thermal conductivity,
and density for different gases.
"""

import chemics as cm
import matplotlib.pyplot as plt

# Parameters
# ----------------------------------------------------------------------------

from params import temp
from params import press

# properties are calculated for each gas item
gas = ['N2', 'H2', 'H2O', 'CO', 'CO2', 'CH4']

# heat capacity [J/(mol K)] for each gas item, same order as `gas` list
# values from Yaw's Handbook at 773.15 K
cp = [31.24, 29.55, 38.32, 31.70, 50.96, 62.13]

# Gas Properties
# ----------------------------------------------------------------------------

mw = []     # store molecular weight of each gas
k = []      # store thermal conductivity of each gas
mu = []     # store viscosity of each gas
rho = []    # store density of each gas

for i in range(len(gas)):
    mw.append(cm.mw(gas[i]))
    mu.append(cm.mu_gas(gas[i], temp))
    rho.append(cm.rhog(mw[i], press, temp))
    if gas[i] == 'CH4':
        k.append(cm.k_gas_organic(gas[i], temp))
    else:
        k.append(cm.k_gas_inorganic(gas[i], temp))


def prandtl(cp, mu, k):
    """
    Prandtl number.

    Parameters
    ----------
    cp : float
        Heat capacity, aka specific heat [J/kgK]
    mu : float
        Gas viscosity [Pa s] or [N s/m²]
    k : float
        Thermal conductivity [W/mK]
    """
    pr = (cp * mu) / k
    return pr


pr = []     # store prandtl number

for i in range(len(cp)):
    c = cp[i] / mw[i] * 1000   # convert J/molK to J/kgK
    m = mu[i] * 1e-6 * 0.1     # convert µP to Ns/m²
    p = prandtl(c, m, k[i])
    pr.append(p)

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

sub = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')

xticks = range(len(gas))
xlabels = [g.translate(sub) for g in gas]

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(8, 6), tight_layout=True)

ax1.bar(xticks, mw, align='center', color='C0')
ax1.set_frame_on(False)
ax1.set_xticks(xticks)
ax1.set_xticklabels(xlabels)
ax1.set_ylabel('MW [g/mol]')

ax2.bar(xticks, mu, align='center', color='C2')
ax2.set_frame_on(False)
ax2.set_xticks(xticks)
ax2.set_xticklabels(xlabels)
ax2.set_ylabel(r'$\mu$ [µP]')

ax3.bar(xticks, rho, align='center', color='C4')
ax3.set_frame_on(False)
ax3.set_xticks(xticks)
ax3.set_xticklabels(xlabels)
ax3.set_ylabel(r'$\rho$ [kg/m³]')

ax4.bar(xticks, k, align='center', color='C5')
ax4.set_frame_on(False)
ax4.set_xticks(xticks)
ax4.set_xticklabels(xlabels)
ax4.set_ylabel('k [W/(m K)]')

ax5.bar(xticks, cp, align='center', color='C1')
ax5.set_frame_on(False)
ax5.set_xticks(xticks)
ax5.set_xticklabels(xlabels)
ax5.set_ylabel('Cp [J/(mol K)]')

ax6.bar(xticks, pr, align='center', color='C9')
ax6.set_frame_on(False)
ax6.set_xticks(xticks)
ax6.set_xticklabels(xlabels)
ax6.set_ylabel('Pr [-]')

plt.show()

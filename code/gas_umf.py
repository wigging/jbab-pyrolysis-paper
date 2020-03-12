"""
Compare minimum fluidization velocity (Umf) of bed material for different
fluidization gases.
"""

import chemics as cm
import matplotlib.pyplot as plt

# Parameters
# ----------------------------------------------------------------------------

# umf is calculated for each gas item
gas = ['N2', 'H2', 'H2O', 'CO', 'CO2', 'CH4']

# pressure [Pa] and temperature [K] of gas in the reactor
p_gas = 101_325
tk_gas = 773.15

# bed particle diameter [m], sphericity [-], and density [kg/m³]
dp = 0.0005
phi = 0.86
rhos = 2500

# void fraction of the bed
ep = 0.46

# Minimum fluidization velocity
# ----------------------------------------------------------------------------

umf_ergun = []
umf_grace = []
umf_wenyu = []

for i in range(len(gas)):
    mw_gas = cm.mw(gas[i])
    mu_gas = cm.mu_gas(gas[i], tk_gas) / 1e7    # convert µP to kg/(ms)
    rho_gas = cm.rhog(mw_gas, p_gas, tk_gas)
    umf_ergun.append(cm.umf_ergun(dp, ep, mu_gas, phi, rho_gas, rhos))
    umf_grace.append(cm.umf_coeff(dp, mu_gas, rho_gas, rhos, coeff='grace'))
    umf_wenyu.append(cm.umf_coeff(dp, mu_gas, rho_gas, rhos, coeff='wenyu'))

# Plot
# ----------------------------------------------------------------------------

sub = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')

xticks = range(len(gas))
xlabels = [g.translate(sub) for g in gas]

fig, ax = plt.subplots(tight_layout=True)
ax.plot(xticks, umf_ergun, marker='.', label='Ergun')
ax.plot(xticks, umf_grace, marker='.', label='Grace')
ax.plot(xticks, umf_wenyu, marker='.', label='WenYu')
ax.grid(color='0.9')
ax.legend(frameon=False, loc='best')
ax.set_frame_on(False)
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels)
ax.set_ylabel('Minimum fluidization velocity, Umf [m/s]')
ax.tick_params(color='0.9')

plt.show()

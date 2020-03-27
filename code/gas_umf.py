"""
Compare minimum fluidization velocity (Umf) of bed material for different
fluidization gases.
"""

import chemics as cm
import matplotlib.pyplot as plt

# Parameters
# ----------------------------------------------------------------------------

from params import dp_bed
from params import ep
from params import phi
from params import press
from params import rhop_bed
from params import temp

# Minimum fluidization velocity
# ----------------------------------------------------------------------------

# umf is calculated for each gas item
gas = ['N2', 'H2', 'H2O', 'CO', 'CO2', 'CH4']

umf_ergun = []
umf_grace = []
umf_wenyu = []

for i in range(len(gas)):
    mw_gas = cm.mw(gas[i])
    mu_gas = cm.mu_gas(gas[i], temp) / 1e7    # convert µP to kg/(ms)
    rho_gas = cm.rhog(mw_gas, press, temp)
    umf_ergun.append(cm.umf_ergun(dp_bed, ep, mu_gas, phi, rho_gas, rhop_bed))
    umf_grace.append(cm.umf_coeff(dp_bed, mu_gas, rho_gas, rhop_bed, coeff='grace'))
    umf_wenyu.append(cm.umf_coeff(dp_bed, mu_gas, rho_gas, rhop_bed, coeff='wenyu'))

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

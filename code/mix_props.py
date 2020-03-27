"""
Compare properties such as molecular weight, viscosity, and density of gas
mixtures.
"""

import chemics as cm
import matplotlib.pyplot as plt

# Parameters
# ----------------------------------------------------------------------------

from params import temp
from params import press

# Mixture properties
# ----------------------------------------------------------------------------

# gas mixtures where each item is a mixture of two gases.
mix_gas = [('H2', 'N2'), ('CO2', 'N2'), ('CO', 'N2'), ('H2O', 'N2')]

# weights for gas mixture where each item is fraction of the two gases [-]
mix_wts = [(0.8, 0.2), (0.7, 0.3), (0.7, 0.3), (0.8, 0.2)]

mw_mix = []     # store molecular weight of each gas mixture
mu_mix = []     # store viscosity of each gas mixture
rho_mix = []    # store density of each gas mixture

for i in range(len(mix_gas)):
    mw1 = cm.mw(mix_gas[i][0])
    mw2 = cm.mw(mix_gas[i][1])
    mu1 = cm.mu_gas(mix_gas[i][0], temp)
    mu2 = cm.mu_gas(mix_gas[i][1], temp)
    wts = mix_wts[i]

    mw_mixture = cm.mw_mix((mw1, mw2), wts)
    mw_mix.append(mw_mixture)

    mu_mixture = cm.mu_herning((mu1, mu2), (mw1, mw2), wts)
    mu_mix.append(mu_mixture)

    rho_mixture = cm.rhog(mw_mixture, press, temp)
    rho_mix.append(rho_mixture)

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

xticks = range(len(mix_gas))
xlabels = ['+'.join(m).translate(sub) for m in mix_gas]

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), tight_layout=True)

ax1.bar(xticks, mw_mix, align='center', color='C0')
for i in range(len(mix_gas)):
    ax1.text(xticks[i] - 0.3, mw_mix[i] + 0.4, f'{mix_wts[i][0]}+{mix_wts[i][1]}', fontsize=9)
ax1.set_frame_on(False)
ax1.set_xticks(xticks)
ax1.set_xticklabels(xlabels)
ax1.set_ylabel('Molecular weight [g/mol]')

ax2.bar(xticks, mu_mix, align='center', color='C4')
for i in range(len(mix_gas)):
    ax2.text(xticks[i] - 0.3, mu_mix[i] + 4, f'{mix_wts[i][0]}+{mix_wts[i][1]}', fontsize=9)
ax2.set_frame_on(False)
ax2.set_xticks(xticks)
ax2.set_xticklabels(xlabels)
ax2.set_ylabel('Viscosity [µP]')

ax3.bar(xticks, rho_mix, align='center', color='C9')
for i in range(len(mix_gas)):
    ax3.text(xticks[i] - 0.3, rho_mix[i] + 0.008, f'{mix_wts[i][0]}+{mix_wts[i][1]}', fontsize=9)
ax3.set_frame_on(False)
ax3.set_xticks(xticks)
ax3.set_xticklabels(xlabels)
ax3.set_ylabel('Density [kg/m³]')

plt.show()

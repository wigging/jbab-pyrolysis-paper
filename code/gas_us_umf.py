"""
Compare minimum fluidization velocity (Umf) of bed material for different
fluidization gases.
"""

import chemics as cm
import matplotlib.pyplot as plt
import numpy as np

# Parameters
# ----------------------------------------------------------------------------

from params import di
from params import dp_bed
from params import ep
from params import phi
from params import press
from params import rhop_bed
from params import temp
from params import q_gas

# Superficial velocity
# ----------------------------------------------------------------------------

ac = (np.pi * di**2) / 4

p_kpa = press / 1000
q_lpm = cm.slm_to_lpm(q_gas, p_kpa, temp)
q_m3s = q_lpm / 60_000
us = q_m3s / ac

# Minimum fluidization velocity
# ----------------------------------------------------------------------------

# umf is calculated for each gas item
gas = ['N2', 'H2', 'H2O', 'CO', 'CO2', 'CH4']

umf_ergun = []
umf_grace = []
umf_rich = []
umf_wenyu = []

us_umf_ergun = []
us_umf_grace = []
us_umf_rich = []
us_umf_wenyu = []

for i in range(len(gas)):
    mw_gas = cm.mw(gas[i])
    mu_gas = cm.mu_gas(gas[i], temp) / 1e7    # convert µP to kg/(ms)
    rho_gas = cm.rhog(mw_gas, press, temp)

    umf_ergun.append(cm.umf_ergun(dp_bed, ep, mu_gas, phi, rho_gas, rhop_bed))
    us_umf_ergun.append(us / umf_ergun[i])

    umf_grace.append(cm.umf_coeff(dp_bed, mu_gas, rho_gas, rhop_bed, coeff='grace'))
    us_umf_grace.append(us / umf_grace[i])

    umf_rich.append(cm.umf_coeff(dp_bed, mu_gas, rho_gas, rhop_bed, coeff='rich'))
    us_umf_rich.append(us / umf_rich[i])

    umf_wenyu.append(cm.umf_coeff(dp_bed, mu_gas, rho_gas, rhop_bed, coeff='wenyu'))
    us_umf_wenyu.append(us / umf_wenyu[i])

# Print
# ----------------------------------------------------------------------------

print(
    f'Parameters\n----------\n'
    f'temp    {temp} K\n'
    f'press   {press:,} Pa\n'
)

print(
    f'Results\n-------\n'
    f'ac      {ac:.4g} m²\n'
    f'us      {us:.4g} m/s\n'
)

print(
    f'Umf [m/s]      {"".join(f"{g:<8}" for g in gas)}\n'
    f'Ergun          {"".join(f"{u:<8.2f}" for u in umf_ergun)}\n'
    f'Grace          {"".join(f"{u:<8.2f}" for u in umf_grace)}\n'
    f'Rich           {"".join(f"{u:<8.2f}" for u in umf_rich)}\n'
    f'WenYu          {"".join(f"{u:<8.2f}" for u in umf_wenyu)}'
)

print(
    f'\n'
    f'Us / Umf [-]   {"".join(f"{g:<8}" for g in gas)}\n'
    f'Ergun          {"".join(f"{u:<8.2f}" for u in us_umf_ergun)}\n'
    f'Grace          {"".join(f"{u:<8.2f}" for u in us_umf_grace)}\n'
    f'Rich           {"".join(f"{u:<8.2f}" for u in us_umf_rich)}\n'
    f'WenYu          {"".join(f"{u:<8.2f}" for u in us_umf_wenyu)}'
)

# Plot
# ----------------------------------------------------------------------------

xticks = np.arange(len(gas))
bar_width = 0.15

sub = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')
xlabels = [g.translate(sub) for g in gas]

fig, ax = plt.subplots(tight_layout=True)
ax.bar(xticks, umf_ergun, color='mediumseagreen', width=bar_width, label='Ergun')
ax.bar(xticks + bar_width, umf_grace, color='mediumblue', width=bar_width, label='Grace')
ax.bar(xticks + bar_width * 2, umf_rich, color='mediumslateblue', width=bar_width, label='Rich')
ax.bar(xticks + bar_width * 3, umf_wenyu, color='dimgrey', width=bar_width, label='WenYu')
ax.legend(frameon=False, loc='best')
ax.set_axisbelow(True)
ax.set_frame_on(False)
ax.set_xticks(xticks + bar_width * 1.5)
ax.set_xticklabels(xlabels)
ax.set_ylabel('Umf [m/s]')
ax.tick_params(bottom=False, left=False)
ax.xaxis.grid(False)
ax.yaxis.grid(True, color='0.9')

fig, ax = plt.subplots(tight_layout=True)
ax.bar(xticks, us_umf_ergun, color='mediumseagreen', width=bar_width, label='Ergun')
ax.bar(xticks + bar_width, us_umf_grace, color='mediumblue', width=bar_width, label='Grace')
ax.bar(xticks + bar_width * 2, us_umf_rich, color='mediumslateblue', width=bar_width, label='Rich')
ax.bar(xticks + bar_width * 3, us_umf_wenyu, color='dimgrey', width=bar_width, label='WenYu')
ax.legend(frameon=False, loc='best')
ax.set_axisbelow(True)
ax.set_frame_on(False)
ax.set_xticks(xticks + bar_width * 1.5)
ax.set_xticklabels(xlabels)
ax.set_ylabel('Us / Umf [-]')
ax.tick_params(bottom=False, left=False)
ax.xaxis.grid(False)
ax.yaxis.grid(True, color='0.9')

plt.show()

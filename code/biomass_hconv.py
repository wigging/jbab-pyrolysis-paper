"""
here
"""

import chemics as cm
import numpy as np

# Parameters
# ----------------------------------------------------------------------------

from params import dp_feed
from params import dp_bed
from params import ep
from params import phi_bed
from params import press
from params import rhop_bed
from params import temp

# Average diameter of biomass particle
# ----------------------------------------------------------------------------

dps = []
wts = []

for x in dp_feed:
    print('d', x['d'], 'mf', x['mf'])
    dps.append(x['d'])
    wts.append(x['mf'])

# biomass particle average Sauter mean diameter [μm]
d_avg = np.average(dps, weights=wts)

# Heat transfer coefficient
# ----------------------------------------------------------------------------

gas = ['N2', 'H2', 'H2O', 'CO', 'CO2', 'CH4']
umf = []
reynolds = []
nusselt = []
hconv = []

for g in gas:
    mw = cm.mw(g)
    rho_gas = cm.rhog(mw, press, temp)
    mu_gas = cm.mu_gas(g, temp) / 1e7    # convert µP to kg/(ms)

    umf_ergun = cm.umf_ergun(dp_bed, ep, mu_gas, phi_bed, rho_gas, rhop_bed)
    umf_grace = cm.umf_coeff(dp_bed, mu_gas, rho_gas, rhop_bed, coeff='grace')
    umf_rich = cm.umf_coeff(dp_bed, mu_gas, rho_gas, rhop_bed, coeff='rich')
    umf_wenyu = cm.umf_coeff(dp_bed, mu_gas, rho_gas, rhop_bed, coeff='wenyu')
    umf_avg = (umf_ergun + umf_grace + umf_rich + umf_wenyu) / 4
    umf.append(umf_avg)

    dp_bio = d_avg / 1e6
    re = (rho_gas * umf_avg * dp_bio) / mu_gas
    reynolds.append(re)

    nu = 2 + (0.9 * re**0.62) * ((dp_bio / dp_bed)**0.2)
    nusselt.append(nu)

    if g == 'CH4':
        k_gas = cm.k_gas_organic(g, temp)
    else:
        k_gas = cm.k_gas_inorganic(g, temp)
    h = (k_gas * nu) / dp_bio
    hconv.append(h)


# Print
# ----------------------------------------------------------------------------

print(
    f'\n{" Parameters ":-^79}\n'
)

print('dp_feed [μm]    mf_feed [%]')
for d, mf in zip(dps, wts):
    print(f'{d:<15} {mf:<15}')

print(
    f'\n{" Results ":-^79}\n\n'
    f'dp feed       {d_avg:.4g} μm (avg.)\n'
    f'dp bed        {dp_bed * 1e6:.1f} μm\n'
)

print(f'{"gas":8} {"Umf":8} {"Re":12} {"Nu":8} {"h":12}')
for i in range(len(gas)):
    print(f'{gas[i]:<8} {umf[i]:<8.2f} {reynolds[i]:<12.2f} {nusselt[i]:<8.2f} {hconv[i]:<12.2f}')

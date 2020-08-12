"""
here
"""

import numpy as np

# Parameters
# ----------------------------------------------------------------------------

from params import dp_feed
from params import dp_bed

# Biomass particle average diameter
# ----------------------------------------------------------------------------

dps = []
wts = []

for x in dp_feed:
    print('d', x['d'], 'mf', x['mf'])
    dps.append(x['d'])
    wts.append(x['mf'])

d_avg = np.average(dps, weights=wts)

# Heat transfer coefficient
# ----------------------------------------------------------------------------

# Print
# ----------------------------------------------------------------------------

print(
    f'\n{" Parameters ":-^79}\n'
)

print('dp_feed [μm]    mf_feed [%]')
for d, mf in zip(dps, wts):
    print(f'{d:<15} {mf:<15}')

print(
    f'\n{" Results ":-^79}\n'
    f'd_avg      {d_avg:.4g} μm\n'
)

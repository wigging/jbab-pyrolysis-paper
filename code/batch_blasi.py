"""
Batch reactor comparing biomass, gas, tar, and char yields using biomass
pyrolysis kinetic reactions from the Di Blasi 1993 and 2001 papers.

References
----------
Colomba Di Blasi. “Analysis of Convection and Secondary Reaction Ef- fects
Within Porous Solid Fuels Undergoing Pyrolysis”. Combustion Science and
Technology, vol. 90, pp. 315–340, 1993.

Colomba Di Blasi and Carmen Branca. “Kinetics of Primary Product Formation
from Wood Pyrolysis”. Industrial & Engineering Chemistry Research, vol. 40,
pp. 5547–5556, 2001.
"""

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# Parameters
# ----------------------------------------------------------------------------

from params import press
from params import temp
from params import y0

# time vector for evaluating kinetic reactions [s]
time = np.linspace(0, 25, num=1000)

# Batch reactor with primary Di Blasi reactions
# ----------------------------------------------------------------------------

gas1 = ct.Solution('blasi.cti')
gas1.TPY = temp, press, y0

# use only primary reactions by disabling the secondary reactions for tar
gas1.set_multiplier(0, 3)    # disable reaction tar => gas
gas1.set_multiplier(0, 4)    # disable reaction tar => char

r1 = ct.IdealGasConstPressureReactor(gas1, energy='off')

sim1 = ct.ReactorNet([r1])
states1 = ct.SolutionArray(gas1, extra=['t'])

for ti in time:
    sim1.advance(ti)
    states1.append(r1.thermo.state, t=ti)

# Batch reactor with primary and secondary Di Blasi reactions
# ----------------------------------------------------------------------------

gas2 = ct.Solution('blasi.cti')
gas2.TPY = temp, press, y0

r2 = ct.IdealGasConstPressureReactor(gas2, energy='off')

sim2 = ct.ReactorNet([r2])
states2 = ct.SolutionArray(gas2, extra=['t'])

for ti in time:
    sim2.advance(ti)
    states2.append(r2.thermo.state, t=ti)

# Batch reactor with primary and secondary Di Blasi reactions (modified)
# ----------------------------------------------------------------------------

gas3 = ct.Solution('blasi.cti')
gas3.TPY = temp, press, y0

# apply factor of 0.2 to reaction tar => gas
gas3.set_multiplier(0.2, 3)

r3 = ct.IdealGasConstPressureReactor(gas3, energy='off')

sim3 = ct.ReactorNet([r3])
states3 = ct.SolutionArray(gas3, extra=['t'])

for ti in time:
    sim3.advance(ti)
    states3.append(r3.thermo.state, t=ti)

# Print
# ----------------------------------------------------------------------------

print(f"""
--- Parameters ---
T   {temp} K
P   {press:,} Pa
""")

print('--- Reactions (index, reaction) ---')
for i, r in enumerate(gas1.reactions()):
    print(i, r)

print('\n--- Final primary yields (mass fraction) ---')
for sp in states1.species_names:
    print(f"{sp:10} {states1(sp).Y[-1][0]:.4f}")

print('\n--- Final primary + seconary yields (mass fraction) ---')
for sp in states2.species_names:
    print(f"{sp:10} {states2(sp).Y[-1][0]:.4f}")

print('\n--- Max tar yield (mass fraction) ---')
print(f"{'tar':10} {max(states1('tar').Y[:, 0]):.4f}   primary")
print(f"{'tar':10} {max(states2('tar').Y[:, 0]):.4f}   primary + secondary")
print(f"{'tar':10} {max(states3('tar').Y[:, 0]):.4f}   primary + secondary (mod)")

print('\n--- Max gas yield (mass fraction) ---')
print(f"{'gas':10} {max(states1('gas').Y[:, 0]):.4f}   primary")
print(f"{'gas':10} {max(states2('gas').Y[:, 0]):.4f}   primary + secondary")
print(f"{'gas':10} {max(states3('gas').Y[:, 0]):.4f}   primary + secondary (mod)")

# Plot
# ----------------------------------------------------------------------------

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 4.8), sharey=True, tight_layout=True)

ax1.plot(states1.t, states1('biomass').Y[:, 0], label='biomass')
ax1.plot(states1.t, states1('gas').Y[:, 0], label='gas')
ax1.plot(states1.t, states1('tar').Y[:, 0], label='tar')
ax1.plot(states1.t, states1('char').Y[:, 0], label='char')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Mass fraction [-]')
ax1.grid(color='0.9')
ax1.set_frame_on(False)
ax1.tick_params(color='0.9')

ax2.plot(states2.t, states2('biomass').Y[:, 0], label='biomass')
ax2.plot(states2.t, states2('gas').Y[:, 0], label='gas')
ax2.plot(states2.t, states2('tar').Y[:, 0], label='tar')
ax2.plot(states2.t, states2('char').Y[:, 0], label='char')
ax2.set_xlabel('Time [s]')
ax2.grid(color='0.9')
ax2.set_frame_on(False)
ax2.tick_params(color='0.9')

ax3.plot(states3.t, states3('biomass').Y[:, 0], label='biomass')
ax3.plot(states3.t, states3('gas').Y[:, 0], label='gas')
ax3.plot(states3.t, states3('tar').Y[:, 0], label='tar')
ax3.plot(states3.t, states3('char').Y[:, 0], label='char')
ax3.set_xlabel('Time [s]')
ax3.grid(color='0.9')
ax3.legend(loc='best', frameon=False)
ax3.set_frame_on(False)
ax3.tick_params(color='0.9')

plt.show()

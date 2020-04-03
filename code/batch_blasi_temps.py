"""
Batch reactor comparing biomass, gas, tar, and char yields for a range of
temperatures using biomass pyrolysis kinetic reactions from the Di Blasi 1993
and 2001 papers.

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

from params import temp_min
from params import temp_max
from params import press
from params import y0

# time vector for evaluating kinetic reactions [s]
time = np.linspace(0, 10, num=1000)

# create range of temperatures [K] in increments of 20
temps = np.arange(temp_min, temp_max + 20, 20)

# Batch reactor with Di Blasi reactions
# ----------------------------------------------------------------------------

# store tar yields at each temperature
# tar1 is for primary reactions only
# tar2 is for primary and secondary reactions
tar1 = []
tar2 = []

# calculate biomass conversion and product yields for each temperature over a
# specified time range
for temp in temps:

    # primary reactions only
    gas1 = ct.Solution('blasi.cti')
    gas1.TPY = temp, press, y0
    gas1.set_multiplier(0, 3)    # disable reaction tar => gas
    gas1.set_multiplier(0, 4)    # disable reaction tar => char

    r1 = ct.IdealGasConstPressureReactor(gas1, energy='off')
    sim1 = ct.ReactorNet([r1])
    states1 = ct.SolutionArray(gas1, extra=['t'])

    # primary and secondary reactions
    gas2 = ct.Solution('blasi.cti')
    gas2.TPY = temp, press, y0

    r2 = ct.IdealGasConstPressureReactor(gas2, energy='off')
    sim2 = ct.ReactorNet([r2])
    states2 = ct.SolutionArray(gas2, extra=['t'])

    for t in time:
        sim1.advance(t)
        states1.append(r1.thermo.state, t=t)

        sim2.advance(t)
        states2.append(r2.thermo.state, t=t)

    tar1.append(states1('tar').Y[:, 0])
    tar2.append(states2('tar').Y[:, 0])

# Print
# ----------------------------------------------------------------------------

print(f"""
--- Parameters ---
temp_min    {temp_min} K
temp_max    {temp_max} K
press       {press:,} Pa
y0          {y0}

--- Calculations ---
temps       {temps} K
""")

# Plot
# ----------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.4, 4.8), sharey=True, tight_layout=True)

for i in range(len(temps)):
    ax1.plot(time, tar1[i])
    ax2.plot(time, tar2[i], label=f'{temps[i]} K')

ax1.grid(color='0.9')
ax1.tick_params(color='0.9')
ax1.set_frame_on(False)
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Tar mass fraction [-]')

ax2.grid(color='0.9')
ax2.legend(loc='best', frameon=False)
ax2.tick_params(color='0.9')
ax2.set_frame_on(False)
ax2.set_xlabel('Time [s]')

plt.show()

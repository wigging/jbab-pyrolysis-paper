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

# Cantera batch reactor example with Blasi biomass pyrolysis kinetics
# ----------------------------------------------------------------------------

# create range of temperatures [K] in increments of 20
temps = np.arange(temp_min, temp_max + 20, 20)

# store biomass conversion and product yields at each temperature
biomass = []
gas = []
tar = []
char = []

# calculate biomass conversion and product yields for each temperature over a
# specified time range
for temp in temps:
    sol = ct.Solution('blasi.cti')
    sol.TPY = temp, press, y0
    r = ct.IdealGasConstPressureReactor(sol, energy='off')

    sim = ct.ReactorNet([r])
    states = ct.SolutionArray(sol, extra=['t'])

    time = np.linspace(0, 10, num=1000)

    for t in time:
        sim.advance(t)
        states.append(r.thermo.state, t=t)

    biomass.append(states('biomass').Y[:, 0])
    gas.append(states('gas').Y[:, 0])
    tar.append(states('tar').Y[:, 0])
    char.append(states('char').Y[:, 0])

# Print
# ----------------------------------------------------------------------------

print(f"""
Parameters
----------
temp_min    {temp_min} K
temp_max    {temp_max} K
press       {press:,} Pa
y0          {y0}

Calculations
------------
temps       {temps} K
""")

# Plot
# ----------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.4, 4.8), tight_layout=True)

for i in range(len(temps)):
    ax1.plot(time, biomass[i], label=f'{temps[i]} K')
    ax2.plot(time, tar[i], label=f'{temps[i]} K')

ax1.grid(color='0.9')
ax1.set_frame_on(False)
ax1.tick_params(color='0.9')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Biomass mass fraction [-]')

ax2.grid(color='0.9')
ax2.legend(loc='best')
ax2.set_frame_on(False)
ax2.tick_params(color='0.9')
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Tar mass fraction [-]')

plt.show()

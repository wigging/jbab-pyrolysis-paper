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

# temperature [K] and pressure [Pa]
tk = 773.15
p = 101_325.0

# initial mass fractions [-]
y = {'biomass': 1, 'gas': 0, 'tar': 0, 'char': 0}

# Cantera batch reactor example with Blasi biomass pyrolysis kinetics
# ----------------------------------------------------------------------------

gas = ct.Solution('blasi.cti')
gas.TPY = tk, p, y
r = ct.IdealGasConstPressureReactor(gas, energy='off')

sim = ct.ReactorNet([r])
states = ct.SolutionArray(gas, extra=['t'])

# time vector to evaluate reaction rates [s]
time = np.linspace(0, 25, num=1000)

for tm in time:
    sim.advance(tm)
    states.append(r.thermo.state, t=tm)

# Print
# ----------------------------------------------------------------------------

print(f"""
Parameters
----------
T {tk} K
P {p:,} Pa
""")

# Plot
# ----------------------------------------------------------------------------

fig, ax = plt.subplots(tight_layout=True)
ax.plot(states.t, states.Y[:, gas.species_index('biomass')], label='biomass')
ax.plot(states.t, states.Y[:, gas.species_index('gas')], label='gas')
ax.plot(states.t, states.Y[:, gas.species_index('tar')], label='tar')
ax.plot(states.t, states.Y[:, gas.species_index('char')], label='char')
ax.grid(color='0.9')
ax.legend(loc='best')
ax.set_frame_on(False)
ax.tick_params(color='0.9')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Mass fraction [-]')

plt.show()

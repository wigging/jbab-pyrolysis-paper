"""
Parameters for the biomass particle, sand particle, and reactor. Units for
each parameter are given in square brackets [].
"""

# Feedstock (biomass particle)
# ----------------------------------------------------------------------------

# particle size distribution of biomass
# 'd' is Sauter mean diameter as μm
# 'mf' is mass fraction as %
dp_feed = [
    {'d': 278, 'mf': 12.1},
    {'d': 344, 'mf': 51.0},
    {'d': 426, 'mf': 34.2},
    {'d': 543, 'mf': 2.7}
]

# Bed (sand particle)
# ----------------------------------------------------------------------------

# sand particle diameter [m]
dp_bed = 0.000453   # 453 μm

# density of a sand particle [kg/m³]
rhop_bed = 2500

# sphericity of a sand particle [-]
phi_bed = 0.94

# Other
# ----------------------------------------------------------------------------

# inner diameter of the bfb reactor [m]
di = 0.05232

# void fraction of the bed [-]
ep = 0.46

# reactor pressure [Pa]
press = 101_325.0

# volumetric flowrate of gas into bfb reactor [SLM]
q_gas = 14

# reactor temperature [K]
temp = 773.15   # 500°C

# minimum and maximum temperatures [K]
temp_min = 753.15   # 480°C
temp_max = 853.15   # 580°C

# initial mass fractions for kinetics [-]
y0 = {'biomass': 1, 'gas': 0, 'tar': 0, 'char': 0}

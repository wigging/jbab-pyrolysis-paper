"""
Parameters for the biomass particle, sand particle, and reactor. Units for
each parameter are given in square brackets [].
"""

# Bed material (sand)
# ----------------------------------------------------------------------------

# bed particle diameter [m]
dp_bed = 0.000453

# density of a bed particle [kg/m³]
rhop_bed = 2500

# sphericity of a bed particle [-]
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

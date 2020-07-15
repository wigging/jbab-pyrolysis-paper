"""
Parameters for the various Python scripts and models.
"""

# inner diameter of the bfb reactor [m]
di = 0.05232

# bed particle diameter [m]
dp_bed = 0.0005

# void fraction of the bed [-]
ep = 0.46

# particle sphericity [-]
phi = 0.86

# reactor pressure [Pa]
press = 101_325.0

# volumetric flowrate of gas into bfb reactor [SLM]
q_gas = 14

# density of a bed particle [kg/m³]
rhop_bed = 2500

# reactor temperature [K]
temp = 773.15   # 500°C

# minimum and maximum temperatures [K]
temp_min = 753.15   # 480°C
temp_max = 853.15   # 580°C

# initial mass fractions for kinetics [-]
y0 = {'biomass': 1, 'gas': 0, 'tar': 0, 'char': 0}

def pyro1(k, kr, rho, cp, r):
    """
    Calculate the pyrolysis number Py I.

    Parameters
    ----------
    k : float
        Thermal conductivity of the biomass particle [W/mK]
    kr : float
        Rate constant [1/s]
    rho : float
        Density of the biomass particle [kg/m³]
    cp : float
        Heat capacity of the biomass particle [J/kgK]
    r : float
        Radius of the biomass particle [m]

    Returns
    -------
    py : float
        Pyrolysis number Py I [-]
    """
    py = k / (kr * rho * cp * (r**2))
    return py

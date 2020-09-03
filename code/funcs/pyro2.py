def pyro2(h, kr, rho, cp, r):
    """
    Calculate the pyrolysis number Py II.

    Parameters
    ----------
    h : float
        Convective heat transfer coefficient [W/m²K]
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
        Pyrolysis number Py II [-]
    """
    py = h / (kr * rho * cp * r)
    return py

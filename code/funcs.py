
def biot(h, r, k):
    """
    Calculate the dimensionless Biot number represented by Bi.

    Parameters
    ----------
    h : float
        Convective heat transfer coefficient [W/m²K]
    r : float
        Radius of the biomass particle [m]
    k : float
        Thermal conductivity of the biomass particle [W/mK]

    Returns
    -------
    Bi : float
        Biot number [-]
    """
    Bi = (h * r) / k
    return Bi


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

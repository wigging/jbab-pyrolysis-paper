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

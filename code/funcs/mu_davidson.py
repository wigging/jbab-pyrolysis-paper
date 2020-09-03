import numpy as np


def mu_davidson(mus, mws, xs):
    """
    Calculate viscosity of a gas mixture.

    Parameters
    ----------
    mus : array_like
        Viscosity of each gas component
    mws : array_like
        Molecular weight of each gas component [g/mol]
    xs : array_like
        Mole fraction of each gas component [-]

    Returns
    -------
    mu : float
        Gas viscosity of the mixture. Units are same as input viscosity.

    Raises
    ------
    here

    Example
    -------
    here

    References
    ----------
    here
    """
    a = 0.375
    e = 2 * np.outer(mws, mws)**0.5 / np.add.outer(mws, mws)
    f = np.sum(np.outer(xs, xs) * e**a / np.outer(mus, mus)**0.5)
    mu_mix = 1 / f
    return mu_mix

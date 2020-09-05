import numpy as np


def mu_davidson(mu, mw, x):
    """
    Calculate viscosity of a gas mixture using the Davidson model [1]_. This
    implementation is based on Equations 26, 28, and 29 in the reference.

    .. math::

       E_{ij} = \\frac{2 \\sqrt{M_i\\, M_j}}{M_i + M_j}

       f = \\sum \\frac{x_i\\, x_j}{\\sqrt{\\mu_i\\, \\mu_j}} E_{ij}^A

       \\mu_{\\text{mix}} = 1 / f

    Parameters
    ----------
    mu : array_like
        Viscosity of each gas component. Units can be µP, cP, μPa·s or some
        other appropriate units for dynamic gas viscosity.
    mw : array_like
        Molecular weight of each gas component [g/mol]
    x : array_like
        Mole fraction of each gas component [-]

    Returns
    -------
    mu_mix : float
        Viscosity of the gas mixture. Units are same as input parameter `mu`.

    Raises
    ------
    ValueError
        If sum of mole fractions does not equal 1.0

    Example
    -------
    Parameters for this example are dynamic gas viscosity in µP, molecular
    weight in g/mol, and mole fraction.

    >>> mu_h2 = 179.75
    ... mu_n2 = 363.87
    ... mw_h2 = 2.016
    ... mw_n2 = 28.014
    ... x_h2 = 0.85
    ... x_n2 = 0.15
    ... mu_wilke([mu_h2, mu_n2], [mw_h2, mw_n2], [x_h2, x_n2])
    206.1662

    References
    ----------
    .. [1] Thomas A. Davidson. A Simple and Accurate Method for Calculating
       Viscosity of Gaseous Mixtures. United States Department of the
       Interior, Report of Investigations 9456, 1993.
    """
    if not np.isclose(np.asarray(x).sum(), 1.0):
        raise ValueError('Sum of mole fractions must be 1.0')

    a = 0.375
    e = 2 * np.outer(mw, mw)**0.5 / np.add.outer(mw, mw)
    f = np.sum(np.outer(x, x) * e**a / np.outer(mu, mu)**0.5)
    mu_mix = 1 / f
    return mu_mix


if __name__ == '__main__':
    # dynamic gas viscosity in µP
    mu_h2 = 179.75
    mu_n2 = 363.87

    # molecular weight in g/mol
    mw_h2 = 2.016
    mw_n2 = 28.014

    # mole fraction
    x_h2 = 0.85
    x_n2 = 0.15

    mu_mix = mu_davidson([mu_h2, mu_n2], [mw_h2, mw_n2], [x_h2, x_n2])
    print(f'mu_mix = {mu_mix:.4f}')

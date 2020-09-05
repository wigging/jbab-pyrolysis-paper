import numpy as np


def mu_wilke(mu, mw, x):
    """
    Calculate viscosity of a gas mixture using approach by Wilke [1]_. This
    implementation is based on Equations 13 and 14 in the reference.

    .. math::

        \\mu_{\\text{mix}} = \\sum_{i=1} \\frac{\\mu_i}{1 + \\frac{1}{x_i} \\sum_{\\substack{j=1 \\j \\neq i}} x_j \\phi_{ij}}

        \\phi_{ij} = \\frac{\\left[1 + (\\mu_i/\\mu_j)^{1/2} (M_j/M_i)^{1/4}\\right]^2}{(4/\\sqrt{2}) \\left[1 + (M_i/M_j)\\right]^{1/2}}

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
    276.4676

    References
    ----------
    .. [1] C.R. Wilke. A Viscosity Equation for Gas Mixtures. The Journal of
       Chemical Physics, vol. 18, no. 4, pp. 517-519, 1950.
    """
    if not np.isclose(np.asarray(x).sum(), 1.0):
        raise ValueError('Sum of mole fractions must be 1.0')

    mi_mj = np.divide.outer(mw, mw)  # Mi / Mj
    num = (1 + (np.divide.outer(mu, mu))**0.5 * (mi_mj.T)**0.25)**2
    den = 4 / np.sqrt(2) * (1 + mi_mj)**0.5
    phi = num / den

    v = x * phi
    vsum = np.sum(v, axis=1) - np.diagonal(v)
    mu_mix = np.sum(mu / (1 + vsum / x))
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
    # x_h2 = 1.0
    # x_n2 = 0.0001

    mu_mix = mu_wilke([mu_h2, mu_n2], [mw_h2, mw_n2], [x_h2, x_n2])
    print(f'mu_mix = {mu_mix:.4f}')

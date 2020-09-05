import numpy as np


def mu_brokaw(mu, mw, x):
    """
    Calculate viscosity of a gas mixture using method by Brokaw [1]_. This
    implementation is for nonpolar gases and is based on Equations 9 and 10 in
    the reference.

    .. math::

       m_{ij} = \\left[ 4 M_i M_j / (M_i + M_j)^2 \\right]^{1/4}

       A_{ij} = m_{ij}\\left( \\frac{M_j}{M_i} \\right)^{1/2} \\left[ 1 + \\frac{\\frac{M_i}{M_j} - \\left(\\frac{M_i}{M_j} \\right)^{0.45} }{ 2 \\left(1 + \\frac{M_i}{M_j} \\right) + \\frac{1 + \\left(\\frac{M_i}{M_j} \\right)^{0.45}}{1 + m_{ij}} m_{ij}} \\right]

       \\mu_{\\text{mix}} = \\sum_{i=1} \\frac{x_i \\sqrt{\\mu_i}}{\\frac{x_i}{\\sqrt{\\mu_i}} + \\sum_{\\substack{j=1 \\ j \\ne i}} \\frac{S_{ij} A_{ij}}{\\sqrt{\\mu_j}} x_j}

    Parameters
    ----------
    mu : array_like
        Viscosity of each gas component. Units can be µP, cP, μPa·s or some
        other appropriate units for dynamic gas viscosity.
    mw : array_like
        Moleculare weight of each gas component [g/mol]
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
    ... mu_brokaw([mu_h2, mu_n2], [mw_h2, mw_n2], [x_h2, x_n2])
    257.9015

    References
    ----------
    .. [1] Richard S. Brokaw. Viscosity of Gas Mixtures. NASA Lewis Research
       Center, NASA technical note NASA-TN-D-4496, 1968.
    """
    if not np.isclose(np.asarray(x).sum(), 1.0):
        raise ValueError('Sum of mole fractions must be 1.0')

    mij = (4 * np.outer(mw, mw) / (np.add.outer(mw, mw)**2))**0.25

    mi_mj = np.divide.outer(mw, mw)  # Mi/Mj
    num = mi_mj - mi_mj**0.45
    den = 2 * (1 + mi_mj) + (1 + mi_mj**0.45) / (1 + mij) * mij
    aij = mij * (mi_mj.T**0.5) * (1 + num / den)

    sij = 1.0
    sqrt_mu = np.sqrt(mu)
    v = sij * aij * x / sqrt_mu
    vsum = np.sum(v, axis=1) - np.diagonal(v)

    mu_mix = np.sum((x * sqrt_mu) / (x / sqrt_mu + vsum))
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

    mu_mix = mu_brokaw([mu_h2, mu_n2], [mw_h2, mw_n2], [x_h2, x_n2])
    print(f'mu_mix = {mu_mix:.4f}')

"""
Calculate minimum fluidization velocity Umf for various gas mixtures.
"""

import chemics as cm

# Get parameters
from params import dp_bed
from params import ep
from params import phi_bed
from params import press
from params import rhop_bed
from params import temp


def run_n2h2():
    """
    Calculate Umf for N2/H2 gas mixtures.
    """
    print('\n--- Umf for N2/H2 gas mixtures ---\n')

    # Mass fractions for N2/H2
    y_n2h2 = [(0.8, 0.2), (0.6, 0.4), (0.4, 0.6), (0.2, 0.8)]

    # Calculate Umf for N2/H2 gas mixtures
    for ys in y_n2h2:

        # Molecular weight for N2 and H2
        mw_n2 = cm.mw('N2')
        mw_h2 = cm.mw('H2')

        # Viscosity for N2 and H2
        mu_n2 = cm.mu_gas('N2', temp) / 1e7
        mu_h2 = cm.mu_gas('H2', temp) / 1e7

        # Molecular weight of mixture
        xs = cm.massfrac_to_molefrac(ys, (mw_n2, mw_h2))
        mw_mixture = cm.mw_mix((mw_n2, mw_h2), xs)

        # Viscosity and density of mixture
        mu_mixture = cm.mu_herning((mu_n2, mu_h2), (mw_n2, mw_h2), xs)
        rho_mixture = cm.rhog(mw_mixture, press, temp)

        # Minimum fluidization velocity of mixture
        umf_mixture = cm.umf_ergun(dp_bed, ep, mu_mixture, phi_bed, rho_mixture, rhop_bed)

        print('ys   ', ys)
        print('xs   ', xs)
        print('umf  ', round(umf_mixture, 4), 'm/s')
        print('')


def run_umf_mix(sp_mix, y):
    """
    Calculate Umf for a given gas mixture `sp_mix` and mass fraction `y`.
    """

    mw1 = cm.mw(sp_mix[0])  # Molecular weight for gas species 1
    mw2 = cm.mw(sp_mix[1])  # Molecular weight for gas species 2

    mu1 = cm.mu_gas(sp_mix[0], temp) / 1e7  # Viscosity for gas species 1
    mu2 = cm.mu_gas(sp_mix[1], temp) / 1e7  # Viscosity for gas species 2

    # Molecular weight of mixture, g/mol
    x = cm.massfrac_to_molefrac(y, (mw1, mw2))
    mw_mixture = cm.mw_mix((mw1, mw2), x)

    # Viscosity and density of mixture
    mu_mixture = cm.mu_herning((mu1, mu2), (mw1, mw2), x)
    rho_mixture = cm.rhog(mw_mixture, press, temp)

    # Minimum fluidization velocity of mixture, m/s
    umf_mixture = cm.umf_ergun(dp_bed, ep, mu_mixture, phi_bed, rho_mixture, rhop_bed)

    print(f'\n--- Umf for {sp_mix} gas mixture ---\n')
    print('y    ', y)
    print('x    ', x)
    print('umf  ', round(umf_mixture, 4), 'm/s')


def run_umf(sp):
    """
    Calculate Umf for a given gas species `sp`.
    """

    mw = cm.mw(sp)                  # Molecular weight, g/mol
    mu = cm.mu_gas(sp, temp) / 1e7  # Gas viscosity, kg/(ms)
    rho = cm.rhog(mw, press, temp)  # Gas density, kg/m³

    # Minimum fluidization velocity, m/s
    umf = cm.umf_ergun(dp_bed, ep, mu, phi_bed, rho, rhop_bed)

    print(f'\n--- Umf for {sp} gas ---\n')
    print('umf  ', round(umf, 4), 'm/s')


def main():
    print('\n--- Parameters ---\n')
    print('dp_bed   ', dp_bed, 'm')
    print('ep       ', ep)
    print('phi_bed  ', phi_bed)
    print('press    ', press, 'Pa')
    print('rhop_bed ', rhop_bed, 'kg/m³')
    print('temp     ', temp, 'K')

    run_n2h2()

    run_umf_mix(('N2', 'CO'), (0.5, 0.5))
    run_umf_mix(('N2', 'CO2'), (0.5, 0.5))

    run_umf('N2')
    run_umf('H2')


if __name__ == '__main__':
    main()

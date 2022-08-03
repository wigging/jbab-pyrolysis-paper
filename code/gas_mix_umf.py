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


def run_n2co():
    """
    Calculate Umf for a N2/CO gas mixture.
    """
    print('\n--- Umf for N2/CO gas mixture ---\n')

    # Mass fraction for N2/CO
    y_n2co = (0.5, 0.5)

    # Molecular weight for N2 and CO
    mw_n2 = cm.mw('N2')
    mw_co = cm.mw('CO')

    # Viscosity for N2 and CO
    mu_n2 = cm.mu_gas('N2', temp) / 1e7
    mu_co = cm.mu_gas('CO', temp) / 1e7

    # Molecular weight of mixture
    xs = cm.massfrac_to_molefrac(y_n2co, (mw_n2, mw_co))
    mw_mixture = cm.mw_mix((mw_n2, mw_co), xs)

    # Viscosity and density of mixture
    mu_mixture = cm.mu_herning((mu_n2, mu_co), (mw_n2, mw_co), xs)
    rho_mixture = cm.rhog(mw_mixture, press, temp)

    # Minimum fluidization velocity of mixture
    umf_mixture = cm.umf_ergun(dp_bed, ep, mu_mixture, phi_bed, rho_mixture, rhop_bed)

    print('ys   ', y_n2co)
    print('xs   ', xs)
    print('umf  ', round(umf_mixture, 4), 'm/s')


def run_n2co2():
    """
    Calculate Umf for a N2/CO2 mixture.
    """
    print('\n--- Umf for N2/CO2 gas mixture ---\n')

    # Mass fraction for N2/CO2
    y_n2co2 = (0.5, 0.5)

    # Molecular weight for N2 and CO2
    mw_n2 = cm.mw('N2')
    mw_co2 = cm.mw('CO2')

    # Viscosity for N2 and CO2
    mu_n2 = cm.mu_gas('N2', temp) / 1e7
    mu_co2 = cm.mu_gas('CO2', temp) / 1e7

    # Molecular weight of mixture
    xs = cm.massfrac_to_molefrac(y_n2co2, (mw_n2, mw_co2))
    mw_mixture = cm.mw_mix((mw_n2, mw_co2), xs)

    # Viscosity and density of mixture
    mu_mixture = cm.mu_herning((mu_n2, mu_co2), (mw_n2, mw_co2), xs)
    rho_mixture = cm.rhog(mw_mixture, press, temp)

    # Minimum fluidization velocity of mixture
    umf_mixture = cm.umf_ergun(dp_bed, ep, mu_mixture, phi_bed, rho_mixture, rhop_bed)

    print('ys   ', y_n2co2)
    print('xs   ', xs)
    print('umf  ', round(umf_mixture, 4), 'm/s')


def run_n2():
    """
    Calculate Umf for N2 gas.
    """
    print('\n--- Umf for N2 gas ---\n')

    # Molecular weight for N2
    mw_n2 = cm.mw('N2')

    # Viscosity for N2
    mu_n2 = cm.mu_gas('N2', temp) / 1e7

    # Density of N2
    rho_n2 = cm.rhog(mw_n2, press, temp)

    # Minimum fluidization velocity
    umf_n2 = cm.umf_ergun(dp_bed, ep, mu_n2, phi_bed, rho_n2, rhop_bed)

    print('umf  ', round(umf_n2, 4), 'm/s')


def run_h2():
    """
    Calculate Umf for H2 gas.
    """
    print('\n--- Umf for H2 gas ---\n')

    # Molecular weight for H2
    mw_h2 = cm.mw('H2')

    # Viscosity for H2
    mu_h2 = cm.mu_gas('H2', temp) / 1e7

    # Density of H2
    rho_h2 = cm.rhog(mw_h2, press, temp)

    # Minimum fluidization velocity
    umf_h2 = cm.umf_ergun(dp_bed, ep, mu_h2, phi_bed, rho_h2, rhop_bed)

    print('umf  ', round(umf_h2, 4), 'm/s')


def main():
    print('\n--- Parameters ---\n')
    print('dp_bed   ', dp_bed, 'm')
    print('ep       ', ep)
    print('phi_bed  ', phi_bed)
    print('press    ', press, 'Pa')
    print('rhop_bed ', rhop_bed, 'kg/m³')
    print('temp     ', temp, 'K')

    run_n2h2()
    run_n2co()
    run_n2co2()
    run_n2()
    run_h2()


if __name__ == '__main__':
    main()

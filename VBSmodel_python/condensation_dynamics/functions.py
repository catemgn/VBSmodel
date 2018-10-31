import scipy.constants as const
import math
import numpy as np
from .default_parameters import dfl


def calculate_mean_molecular_speed(molecular_mass, temperature=dfl['values'][
    'temperature']):
    """
    Calculates the mean molecular speed of a organic gas molecule [m/s].(It is
    related to the collision speed between gas molecule and aerosol particle).

    :param molecular_mass:
         mass of the organic gas molecule [g/mol].
    :type molecular_mass: float

    :param temperature:
         ambient temperature considered [K].
    :type temperature: float

    :return:
        mean molecular speed [m/s].
    :rtype: float
    """

    return math.sqrt((8 * dfl['conversions'][
        'kg2g'] * const.R * temperature) / (const.pi * molecular_mass))


def calculate_molecular_diffusion_constant(molecular_mass=dfl['values'][
    'gas_molar_weight'], gas_density=dfl['values']['gas_density'],
                                           temperature=dfl['values'][
                                               'temperature2'],
                                           pressure=dfl['values'][
                                               'pressure']):
    """
    Calculates the diffusion constant of the organic gas molecule [m2/s].

    :param molecular_mass:
        mass of the gas molecule [g/mol].
    :type molecular_mass: float

    :param: gas_density:
        typical organic gas density [kg/m3].
    :type gas_density: float

    :param temperature:
        ambient temperature considered [K].
    :type temperature: float

     :param pressure:
        ambient pressure considered [Pa].
    :type pressure: float

    :return:
        gas molecular diffusion constant D.
    :rtype: float
    """

    # Scaling to a reference compound values: (#TODO to put in default later?).
    d_ref = 6.6e-6  # reference diffusion constant [m2/s].
    t_ref = 300  # [K] reference temperature.
    p_ref = 101325  # [Pa] reference pressure.
    m_ref = 200  # [g / mol] reference molecular weight.
    rho_ref = 1200  # [kg/m3] reference molecular density.

    # Other parameters needed
    m_air = dfl['values']['air_mass']
    rho_air = dfl['values']['air_density']
    m = dfl['functions']['molecular_diffusion_constant']['m']
    delta = dfl['functions']['molecular_diffusion_constant']['delta']

    # renaming the input variable for simplicity.
    m_i = molecular_mass
    rho_i = gas_density

    return d_ref * (temperature / t_ref) ** m * (pressure / p_ref) * math.sqrt(
        m_ref * (m_i + m_air) / (m_i * (
                m_ref + m_air))) * ((m_ref / rho_ref + m_air / rho_air) / (
            m_i / rho_i + m_air / rho_air)) ** delta


def calculate_particle_critical_diameter(mean_molecular_speed,
                                         diffusion_constant):
    """
    Calculates the critical diameter, i.e. the crossover between kinetic and
    continuum deposition speeds as defined in Donahue 2017 [nm].

    :param mean_molecular_speed:
        mean molecular speed of a organic gas molecule [m/s].
    :type mean_molecular_speed: float

    :param diffusion_constant:
        diffusion constant of the organic gas molecule [m2/s].
     :type diffusion_constant: float

    :return:
        particle critical diameter [nm].
    :rtype: float
    """

    return 8 * diffusion_constant / (mean_molecular_speed * dfl['conversions'][
        'nm2m'])


def calculate_particle_volume(particle_diameter):
    """
    Calculates the particle spherical equivalent volume [m3].

    :param particle_diameter:
        particle diameter [nm].
    :type particle_diameter: float

    :return:
        particle volume [m3].
    :rtype: float
    """

    return math.pi / 6 * (particle_diameter * dfl['conversions'][
        'nm2m']) ** 3


def calculate_particle_mass(particle_volume, gas_density=dfl['values'][
    'gas_density']):
    """
    Calculates the particle mass [g/mol].

    :param particle_volume:
        volume of the article [m3].
    :type particle_volume; float

    :param gas_density:
        density of the gas molecules [kg/m3].
    :type gas_density; float

    :return:
        particle mass [g/mol].
    :rtype: float
    """

    return particle_volume * gas_density / dfl['conversions']['g/mol2kg']


def calculate_reduced_mass(particle_mass, molecular_mass):
    """
    Calculates the reduced mass of the collision (vapor and particle) [g/mol].

    :param particle_mass:
         mass of the particle [g/mol].
    :type particle_mass: float.

    :param molecular_mass:
        mass of the organic gas molecule [g/mol].
    :type molecular_mass: float

    :return:
        reduced mass [g/mol].
    :rtype: float
    """

    return particle_mass * molecular_mass / (particle_mass + molecular_mass)


# TODO this is a repetition of function calculate_mean of molecular_speed.
# TODO Find a more elegant way to make them in one single function.
def calculate_center_mass_speed_1(reduced_mass, temperature=dfl[
    'values']['temperature']):
    """
    Calculates the center of mass speed of the binary system gas+ particle[m/s].

    :param reduced_mass:
         reduced mass of the system gas+particle [g/mol].
    :type reduced_mass: float

    :param temperature:
         ambient temperature considered [K].
    :type temperature: float

    :return:
        center of mass speed of the system gas+particle [m/s].
    :rtype: float
    """

    return (8 * dfl['conversions']['kg2g'] * const.R * temperature / (
            const.pi * reduced_mass)) ** 0.5 / 4


def calculate_speed_enhancement_factor(particle_mass, molecular_mass):
    """
     Calculates the speed enhancement factor of the collision [-].

     :param particle_mass:
          mass of the particle [g/mol].
     :type particle_mass: float.

     :param molecular_mass:
         mass of the organic gas molecule [g/mol].
     :type molecular_mass: float

     :return:
         speed enhancement factor [-].
     :rtype: float
     """

    return ((particle_mass + molecular_mass) / particle_mass) ** 0.5


def calculate_center_mass_speed_2(mean_molecular_speed,
                                  speed_enhancement_factor):
    """
    Calculates the center of mass speed of the system gas+particle (alternative way) [m/s].
    :param mean_molecular_speed:
        mean speed of the gas molecule [m/s].
    :type mean_molecular_speed: float

    :param speed_enhancement_factor:
        speed enhancement factor of the collision [-].
    :type speed_enhancement_factor: float

    :return:
        center of mass speed of the system gas+particle [m/s].
    :rtype:
    """
    return mean_molecular_speed * speed_enhancement_factor


def calculate_molecular_size_enhancement(particle_diameter, molecular_diameter):
    """
    Calculates the molecular enhancement, to relate the circular cross
    section to the spherical surface area [-].

    :param particle_diameter:
        diameter of the particle [nm].
    :type particle_diameter: float

    :param molecular_diameter:
        diameter of the gas molecule [nm].
    :type molecular_diameter: float

    :return:
        molecular size enhancement [-].
    :rtype: float
    """

    return (particle_diameter + molecular_diameter) ** 2 / (
            particle_diameter ** 2)


def calculate_knudsen_number(particle_diameter, particle_critical_diameter):
    """
    calculates a modified Knudsen Number based on the critical particle
    diameter [-].

    :param particle_diameter:
        diameter of the particle [nm].
    :type particle_diameter: float

    :param particle_critical_diameter:
        critical diameter of the particle [nm].
    :type particle_critical_diameter: float

    :return:
        Knudsen Number [-].
    :rtype: float
    """

    # TODO why is 3/4? not in the paper...
    return 3 / 4 * particle_critical_diameter / particle_diameter


# TODO revise formulas inside the transition regime function (bsxfun)
def calculate_transition_regime_correction(accommodation_coefficient,
                                           knudsen_number):
    """
    calculates the transition regime correction factor based on Fuchs and
    Sutugun.

    :param accommodation_coefficient:
        accommodation coefficient of the gas molecule on the particle [-].
    :type accommodation_coefficient: float

    :param knudsen_number:
        knudsen number [-].
    :type knudsen_number: float

    :return:
        transition regime correction from kinetic to continuum [-].
    :rtype: float
    """

    # TODO why is defined like that check the paper...?
    # calculate the Knudsen Number:
    # km = (3 / 4) * particle_diameter / critical_particle_diameter
    km = knudsen_number
    return km * (1 + km) / (km ** 2 + km + 0.75 * accommodation_coefficient +
                            0.283 * km * accommodation_coefficient)


def calculate_full_deposition_speed(accommodation_coefficient,
                                    molecular_size_enhancement,
                                    transition_correction, cm_speed):
    """
    :param accommodation_coefficient:
        accommodation coefficient of the gas molecule on the particle [-].
    :type accommodation_coefficient: float

    :param molecular_size_enhancement:
        molecular size enhancement [-].
    :type molecular_size_enhancement: float

    :param transition_correction:
        transition regime correction from kinetic to continuum [-].
    :type transition_correction: float

    :param cm_speed:
        center of mass speed (gas + particle) [m/s].
    :type cm_speed: float

    :return:
        full corrected deposition speed [m/s].
        :rtype: float
    """

    return accommodation_coefficient * molecular_size_enhancement * cm_speed * \
           transition_correction


def calculate_kinetic_deposition_speed(mean_molecular_speed,
                                       accommodation_coefficient):
    """
    Calculates the deposition speed at the kinetic limit [m/s].

    :param mean_molecular_speed:
        mean molecular speed of a organic gas molecule [m/s].
    :type mean_molecular_speed: float
    :param accommodation_coefficient:
        accommodation coefficient of the gas molecule on the particle [-].
    :type accommodation_coefficient: float

    :return:
     kinetic deposition speed [m/s].
    :rtype: float
    """

    return mean_molecular_speed * accommodation_coefficient / 4


def calculate_continuum_deposition_speed(particle_diameter,
                                         diffusion_constant):
    """
    Calculates the deposition speed at the continuum limit [m/s].

    :param diffusion_constant:
        diffusion constant of the organic gas molecule [m2/s].
     :type diffusion_constant: float

    :param particle_diameter:
        particle diameter [nm].
    :type particle_diameter: float

    :return:
     continuum deposition speed [m/s].
    :rtype: float
    """

    return 2 * diffusion_constant / (
            particle_diameter * dfl['conversions']['nm2m'])


def calculate_deposition_speed_normalised_ratio(full_deposition_speed,
                                                kinetic_deposition_speed):
    """
    calculates the normalized ratio of the full deposition speed to the
    kinetic deposition speed [-].

    :param full_deposition_speed:
        full corrected deposition speed [m/s].
    :type full_deposition_speed: float

    :param kinetic_deposition_speed:
        deposition speed in the kinetic limit [m/s].
    :type kinetic_deposition_speed: float

    :return:
        normalised ratio of speed deposition [-].
    :rtype: float
    """

    return full_deposition_speed / kinetic_deposition_speed


def calculate_diameter_growth_rate(full_deposition_speed, gas_density):
    """
    Calculates the rate of diameter change per unit condensation driving
    force (i.e. a driving force of 1ug/m3) [nm/hr].

    :param full_deposition_speed:
        full corrected deposition speed [m/s].
    :type full_deposition_speed: float

    :param gas_density:
        density of the gas molecules [kg/m3].
    :type gas_density; float

    :return:
        diameter growth rate [nm/hr].
    :rtype: float
    """

    return 2 * full_deposition_speed / gas_density / dfl['conversions'][
        'sec2hr']


def calculate_kelvin_term(particle_diameter, kelvin_diameter):
    """
    Calculates the Kelvin Term [-]. #TODO improve function documentation.

    :param particle_diameter:
        particle diameter [nm].
    :type particle_diameter: float

    :param kelvin_diameter:
        Kelvin diameter [nm].
    :type kelvin_diameter: float

    :return:
        Kelvin Term [-].
    :rtype: float
    """

    return 10 ** (kelvin_diameter / particle_diameter)


def calculate_evaporation_timescale(particle_diameter, kelvin_term,
                                    full_deposition_speed,
                                    saturation_concentration, gas_density):
    """
    Calculates the time for evaporation of constituent from a particle
    at the low mass fraction limit (i.e. when the diameter does not change as
    the constituent evaporates [hr].

    :param particle_diameter:
        particle diameter [nm].
    :type particle_diameter: float

    :param full_deposition_speed:
        full corrected deposition speed [m/s].
    :type full_deposition_speed: float

    :param kelvin_term:
        kelvin_term [-].
    :type kelvin_term: float

    :param: gas_density:
        typical organic gas density [kg/m3].
    :type gas_density: float

    :param saturation_concentration:
        saturation concentration (volatility) [ug/m3].
    :type saturation_concentration: float

    :return:
        evaporation timescale [hr].
    :rtype: float
    """

    return dfl['conversions']['sec2hr'] * (
                gas_density * particle_diameter * dfl['conversions'][
            'nm2m']) / (6 * full_deposition_speed * kelvin_term * saturation_concentration /
                   dfl['conversions']['kg2ug'])


def calculate_particle_surface(particle_diameter):
    """
    Calculates the surface area of a particle for a given diameter [nm2].

    :param particle_diameter:
        particle diameter [nm].
    :type particle_diameter: float

    :return:
        particle surface area [nm2].
    :rtype:float
    """

    return math.pi * particle_diameter ** 2


# TODO remember that for coll freq the dep speed is calculated with alpha =1!
def calculate_collision_frequency(particle_surface, full_deposition_speed,
                                  molecular_mass,
                                  vapour_concentration=dfl['values']['Cv']):
    """
    Calculates the frequency with which molecules collide with a particle [1/s].

    :param particle_surface:
        surface area of the particle [nm2].
    :type particle_surface: float

    :param full_deposition_speed:
        full corrected deposition speed [m/s].
    :type full_deposition_speed: float

    :param molecular_mass:
         mass of the organic gas molecule [g/mol].
    :type molecular_mass: float

    :param vapour_concentration:
        Concentration of the species in the gas phase [ug/m3].
    :type vapour_concentration: float

    :return:
        molecule-particle collision frequency [1/s].
    :rtype: float
    """
    #TODO check conversions...
    return particle_surface * (dfl['conversions'][
        'nm2m'] ** 2) * full_deposition_speed * (vapour_concentration /dfl[
        'conversions']['kg2ug'])/(molecular_mass*dfl['conversions'][
        'g/mol2kg'])

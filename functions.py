import numpy as np
import scipy.constants as scicon
import math
from default_parameters import dfl
import utils as utl


PRODUCTION_THRESHOLD = -20  # Time threshold for stop net SOA production [min].
NMtoM = 1e-9
KGM3toUGM3 = 1e9
KGtoUG = 1e9
CM3toM3 = 1e-6
SECtoMIN = 60
HOURtoSEC = 3600
GtoKG = 1e-3
AMUtoKG = 1.660538921e-27
UGtoAMU = 1e-9 / AMUtoKG


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

    return ((8 * dfl['conversions'][
        'kg2g'] * scicon.R * temperature) / (scicon.pi * molecular_mass)) ** 0.5


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

    return d_ref * (temperature / t_ref) ** m * (pressure / p_ref) * ((
            m_ref * (m_i + m_air) / (m_i * (
            m_ref + m_air)) )** 0.5) * ((m_ref / rho_ref + m_air / rho_air) / (
            m_i / rho_i + m_air / rho_air)) ** delta


def calculate_particle_critical_diameter(mean_molecular_speed,
                                         diffusion_constant,
                                         accommodation_coefficient):
    """
    Calculates the critical diameter, i.e. the crossover between kinetic and
    continuum deposition speeds as defined in Donahue 2017 [nm].

    :param mean_molecular_speed:
        mean molecular speed of a organic gas molecule [m/s].
    :type mean_molecular_speed: float

    :param diffusion_constant:
        diffusion constant of the organic gas molecule [m2/s].
    :type diffusion_constant: float

    :param accommodation_coefficient:
        accommodation coefficient of the gas molecule on the particle [-].
    :type accommodation_coefficient: float

    :return:
        particle critical diameter [nm].
    :rtype: float
    """

    return 8 * diffusion_constant / (
            accommodation_coefficient * mean_molecular_speed * dfl[
        'conversions']['nm2m'])


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

    mi_dot_mp = utl.__matrix_multiply(molecular_mass, particle_mass)
    a = utl.__extend_to_matrix(molecular_mass, len(particle_mass))
    mi_plus_mp = utl.__column_by_vector(a, particle_mass, '+')

    return mi_dot_mp / mi_plus_mp


# TODO this is a repetition of function calculate_mean of molecular_speed.
# TODO Find a more elegant way to make them in one single function (utility?).
def calculate_center_mass_speed(reduced_mass, temperature=dfl[
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

    return (8 * dfl['conversions']['kg2g'] * scicon.R * temperature / (
            scicon.pi * reduced_mass)) ** 0.5


def calculate_molecular_size_enhancement(particle_diameter, molecular_diameter):
    """
    Calculates the molecular enhancement, to relate the circular cross
    section to the spherical surface area [-].

    :param particle_diameter:
        diameter of the particle [nm].
    :type particle_diameter: float | np.array

    :param molecular_diameter:
        diameter of the gas molecule [nm].
    :type molecular_diameter: float | np.array

    :return:
        molecular size enhancement [-].
    :rtype: float or np.array
    """

    a = utl.__extend_to_matrix(molecular_diameter, len(particle_diameter))
    di_plus_dp = utl.__column_by_vector(a, particle_diameter, '+')

    return utl.__column_by_vector(di_plus_dp**2, particle_diameter**2, '/')


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

    a = utl.__extend_to_matrix(particle_critical_diameter,
                               len(particle_diameter))

    return utl.__column_by_vector(a, particle_diameter, '/')




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

    alpha = accommodation_coefficient
    km = knudsen_number
    up = (km + 1) * (0.75 * alpha * km + 1) / (
            km * (0.75 * alpha * km + 1 + 0.283 * alpha) + 1)

    return km / (km + 1) * up


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

    return accommodation_coefficient * molecular_size_enhancement * (
            cm_speed / 4) * transition_correction


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
    const = 2/NMtoM
    a = utl.__extend_to_matrix(diffusion_constant, len(particle_diameter))

    return const * utl.__column_by_vector(a, particle_diameter, '/')


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

    return utl.__row_by_vector(full_deposition_speed,
                               kinetic_deposition_speed, '/')


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


def calculate_kelvin_terms_Kp(particle_diameter, kelvin_diameter):
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
    conv_factor = 1 / HOURtoSEC * NMtoM * KGtoUG
    num = (gas_density * particle_diameter)/6
    a = utl.__matrix_multiply(saturation_concentration, kelvin_term)
    den = np.multiply(full_deposition_speed, a)

    return conv_factor*num/den


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
def calculate_collision_frequency_on_particles(particle_surface,
                                               full_deposition_speed,
                                               molecular_mass,
                                               vapour_concentration=
                                               dfl['values']['Cv']):
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
    conv_factor = NMtoM ** 2 / KGtoUG / AMUtoKG
    # TODO check conversions...
    # return conv_factor * particle_surface * full_deposition_speed * (
    # vapour_concentration / molecular_mass)

    a = full_deposition_speed * particle_surface[:, None].T
    b = vapour_concentration / molecular_mass
    return conv_factor * a * b[:, None]


def calculate_concentrations_Cs_i(concentrations_Cs_ip):
    """
    Calculates the suspended concentration for each species i Cs_ip [ug/m3].

    :param concentrations_Cs_ip:
        suspended concentrations for each species i and population p [ug/m3].
    :type concentrations_Cs_ip: numpy.array (2D).

    :return:
        suspended concentrations for each species i [ug/m3].
    :rtype: numpy.array (1D).
    """
    return concentrations_Cs_ip.sum(axis=1)  # sum over rows.


def calculate_concentrations_Cs_p(concentrations_Cs_ip):
    """
    Calculates the suspended concentration for each population p Cs_p ug/m3].

    :param concentrations_Cs_ip:
        suspended concentrations for each species i and population p [ug/m3].
    :type concentrations_Cs_ip: numpy.array (2D).

    :return:
        suspended concentrations for each population p [ug/m3].
    :rtype: numpy.array (1D).
    """
    return concentrations_Cs_ip.sum(axis=0)


def calculate_concentration_Cs_OA(concentrations_Cs_i):
    """
    Calculates the total suspended concentration Cs_OA [ug/m3].

    :param concentrations_Cs_i:
        suspended concentrations for each species i [ug/m3].
    :type concentrations_Cs_i: numpy.array (1D).

    :return:
       total suspended concentration [ug/m3].
    :rtype: float.
    """
    return np.sum(concentrations_Cs_i)


def calculate_concentrations_Cs_seed_p(seed_density_rho_seed_p,
                                       seed_diameters_d_seed_p,
                                       number_concetrations_Ns_p):
    """
    Calculates the seed concentration for each population p [ug/m3].

    :param seed_density_rho_seed_p:
        seed density for each population p [kg/m3].
    :type seed_density_rho_seed_p: numpy.array (1D).

    :param number_concentrations_Ns_p:
        suspended number concentrations for each population p [1/m3].
    :type number_concentrations_Ns_p: numpy.array (1D).

    :param seed_diameters_d_seed_p:
        seed diameter for each population p.
    :type seed_diameters_d_seed_p: numpy.array (1D).

    :return: seed concentrations C_seed_p [ug/m3].
    :rtype: numpy.array (1D).
    """
    conv_fact = scicon.pi / 6 * NMtoM**3 * KGM3toUGM3 /CM3toM3
    return conv_fact * (seed_diameters_d_seed_p) ** 3 * (
            seed_density_rho_seed_p) * (
                   number_concetrations_Ns_p)


def calculate_concentration_Cs_seed(concentrations_Cs_seed_p):
    """
    Calculates the total suspended seed concentration Cs_seed [ug/m3].

    :param concentrations_Cs_seed_p:
        suspended seed concentrations for each population p [ug/m3].
    :type concentrations_Cs_seed_p: numpy.array (1D).

    :return:
        total suspended seed concentration [ug/m3].
    :rtype: float.
    """
    return np.sum(concentrations_Cs_seed_p)


def calculate_number_concentration_Ns(number_concentrations_Ns_p):
    """
    Calculates the total suspended number concentration Ns [ 1/m3].

    :param number_concentrations_Ns_p:
        suspended number concentrations for each population p [1/m3].
    :type number_concentrations_Ns_p: numpy.array (1D).

    :return:
        total suspended number concentration [1/m3].
    :rtype: float.
    """
    return np.sum(number_concentrations_Ns_p)


def set_precursor_loss_Lprec(production_time, precursor_loss_Lprec):
    """
    Stops the constant production of SOA given a time threshold.
    (NB: the SOA production is linked to the precursor chemical loss).

    :param production_time:
        interval of constant SOA production [min].
    :type production_time: numpy.array (1D).

     :param precursor_loss_Lprec:
        precursor loss rate  [ug/(m3*min)].
    :type precursor_loss_Lprec: float.

    :return: precursor loss rate [ug/(m3*min)].
    :rtype: float.
    """

    if production_time > PRODUCTION_THRESHOLD:  # TODO CHECK IF CORRECT WHEN
        # TEST
        precursor_loss_Lprec = 0

    return precursor_loss_Lprec


def calculate_net_vapour_productions_Pv_i(fraction_yields_yi,
                                          precursor_loss_Lprec):
    """
    Calculates the net vapour production for each species i [-].
    :param fraction_yields_yi:
        fractions of yield for each species i [-].
    :type fraction_yields_yi: numpy.array (1D).

    :param precursor_loss_Lprec:
        precursor loss rate [ug/(m3*min)].

    :return: net_vapour_productions_Pv_i [ug/m3].
    :rtype: numpy.array (1D).
    """

    return precursor_loss_Lprec * fraction_yields_yi


def calculate_activities_as_ip(concentrations_Cs_ip, concentrations_Cs_p):
    """
    Calculates the activity of each species i and population p [-].

    :param concentrations_Cs_ip:
        suspended concentrations for each species i and population p [ug/m3].
    :type concentrations_Cs_ip: numpy.array (2D).

    :param concentrations_Cs_p:
        suspended concentrations for each population p [ug/m3].
    :type concentrations_Cs_p: numpy.array (1D).

    :return: activity of each species i and population p [-].
    :rtype: numpy.array (2D).
    """
    is_concentration = concentrations_Cs_p != 0  # check if there are 0 values.
    eps =2.2204e-16
    return is_concentration * utl.__column_by_vector(concentrations_Cs_ip,
                                                     concentrations_Cs_p +eps,
                                                     '/')


def calculate_condensation_driving_forces_Fvs_ip(activities_as_ip,
                                                 saturation_concentrations_Co_i,
                                                 concentrations_Cv_i,
                                                 kelvin_terms_Kp):
    """
    Calculates the condensation driving force from vapour (v) to suspended (s)
    phases for each species i and population p [ug/m3].

    :param: saturation_concentrations_Co_i:
           saturation concentration for each species i (VBS set) [ug/m3].
    :type: saturation_concentrations_Co_i: numpy.array (1D).

    :param: activities_as_ip:
        activity of each species i and population p [-].
    :type activities_as_ip: numpy.array (2D).

    :param: concentrations_Cv_i:
        vapour concentrations for each species i[ug/m3].
    :type: numpy.array (1D).

    :param: kelvin_terms_Kp:
        kelvin term for each population p [-].
    :type: numpy.array (1D).

    :return: condensation_driving_forces for each species i and population p
    [ug/m3].
    :rtype: numpy.array (2D).
    """

    a = utl.__row_by_vector(activities_as_ip, saturation_concentrations_Co_i, '*')
    b = utl.__column_by_vector(a, kelvin_terms_Kp, '*')
    return utl.__row_by_vector(-b, concentrations_Cv_i, '+')


def calculate_particle_volumes_vs_p(number_concentrations_Ns_p,
                                    seed_density_rho_seed_p, organics_density,
                                    concentrations_Cs_p,
                                    concentrations_Cs_seed_p):
    """
    Calculates the particle volume for each population p [m3].

    :param number_concentrations_Ns_p:
        suspended number concentrations for each population p [1/m3].
    :type number_concentrations_Ns_p: numpy.array (1D).

    :param seed_density_rho_seed_p:
        density of seed particles for each population pm[ug/m3].
    :type seed_density_rho_seed_p: np.array

     :param organics_density:
        density of organic particles [ug/m3].
    :type organics_density: float.

    :param concentrations_Cs_p:
        suspended concentrations for each population p [ug/m3].
    :type concentrations_Cs_p: numpy.array (1D).

    :param concentrations_Cs_seed_p:
        suspended seed concentrations for each population p [ug/m3].
    :type concentrations_Cs_seed_p: numpy.array (1D).

    :return: particle volume for each population p [m3].
    :rtype: numpy.array (1D).

    """
    conv_factor = CM3toM3/KGM3toUGM3
    return conv_factor*((concentrations_Cs_seed_p / seed_density_rho_seed_p)
            + (concentrations_Cs_p /
               organics_density)) / number_concentrations_Ns_p


def calculate_particle_diameters_ds_p(particle_volumes_vs_p):
    """
    Calculates the particle diameter for each population p [m].

    :param particle_volumes_vs_p:
        particle volume for each population p [m3].
    :type particle_volumes_vs_p: numpy.array (1D).

    :return: particle diameter for each population p [m3].
    :rtype: numpy.array (1D)
    """

    return np.cbrt((6 / scicon.pi * particle_volumes_vs_p))/NMtoM


def calculate_collision_frequencies_nus_ip(full_deposition_speed,
                                           particle_diameters_ds_p,
                                           number_concentrations_Ns_p):
    """
    Calculates the collision frequency of vapour i over particle in
    population p [1/min].

    :param full_deposition_speed:
        full corrected deposition speed [m/s].
    :type full_deposition_speed: numpy.array (1D).

    :param particle_diameters_ds_p:
        particle diameter for each population p [m3].
    :type particle_diameters_ds_p: numpy.array (1D).

     :param number_concentrations_Ns_p:
        suspended number concentrations for each population p [1/m3].
    :type number_concentrations_Ns_p: numpy.array (1D).

    :return: collision frequency between vapour i and particle population p.
    :rtype: numpy.array (2D).
    """
    conv_fact = scicon.pi*NMtoM**2*SECtoMIN/CM3toM3
    a =np.square(particle_diameters_ds_p) * number_concentrations_Ns_p

    return conv_fact * utl.__column_by_vector(full_deposition_speed, a, '*')


def calculate_condensation_sinks_ks_ip(collision_frequency_nus_ip,
                                       effective_accommodation_coefficient=1):
    """
    Calculates the condensation sinks for vapour of species i on
    suspended population p [1/min].

    :param collision_frequency_nus_ip:
        collision frequency between vapour i and particle population p [1/min].
    :type collision_frequency_nus_ip: numpy.array (1D).

    :param: efficient_accommodation_coeff
    :type float

    :return: condensation sinks for vapour of species i on population p [1/min].
    :rtype: numpy.array (2D).
    """

    return effective_accommodation_coefficient *collision_frequency_nus_ip


def calculate_condensation_fluxes_phi_vs_ip(condensation_sinks_ks_ip,
                                            condensation_driving_forces_Fvs_ip):
    """
    Calculates the condensations fluxes from vapour (v) to suspended (s)
    phase for each species i and population p [ug/m3].

    :param: condensation_sinks_ks_ip:
        condensation sinks for vapour of species i on population p [1/min].
    :type: numpy.array (1D).

    :param condensation_driving_forces_Fvs_ip:
        condensation driving force for each species i and population p [ug/m3].
    :type: numpy.array (2D).

    :return: condensation fluxes for each species i and population p [ug/m3].
    :rtype: numpy.array (2D).
    """

    return np.multiply(condensation_sinks_ks_ip,
                       condensation_driving_forces_Fvs_ip)


def calculate_vapours_change_rates_dCv_i(condensation_fluxes_phi_vs_ip,
                                         net_vapour_productions_Pv_i):
    """
    Calculates the vapour concentration change rate for each species i
    dCv_i/dt [ug/(m3*min)].

    :param condensation_fluxes_phi_vs_ip:
        fluxes for each species i and population p [ug/m3].
    :type: numpy.array (2D).

    :param net_vapour_productions_Pv_i:
        net vapour production for each species i [ug/(m3*min)].
    :type: numpy.array (1D).

    :return: vapour concentration change rate for each species i.
    :rtype: numpy.array (1D).
    """

    return net_vapour_productions_Pv_i - condensation_fluxes_phi_vs_ip.sum(
        axis=1)


def calculate_suspended_particles_change_rates_dCs_ip(
        condensation_fluxes_phi_vs_ip,
        net_particle_productions_Ps_ip):
    """
    Calculates the suspended concentration change rate for each species i and population p
    dCs_ip/dt [ug/(m3*min)].

    :param condensation_fluxes_phi_vs_ip:
        fluxes for each species i and population p [ug/m3].
    :type: numpy.array (2D).

    :param net_particle_productions_Ps_ip:
        net suspended particle production for each species i and population p [
        ug/m3].
    :type: numpy.array (2D).

    :return: suspended particle concentration change rate for each species i
    and population p.
    :rtype: numpy.array (2D).
    """

    return net_particle_productions_Ps_ip + condensation_fluxes_phi_vs_ip







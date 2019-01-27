import schedula as sh
import functions as fun


def model():
    """
    Returns the model dispatcher.
    WARNING! make sure that the order of the inputs is exactly the same of the
    order of inputs in the corresponding function.

    :return:
        workflow map of the model.
    :rtype: schedula.Dispatcher
    """

    dsp = sh.Dispatcher(
        name='condensation_map',
        description='Provides the Condensation Dynamics of the VBS model.'
    )

    dsp.add_function(
        function=fun.calculate_mean_molecular_speed,
        inputs=['molecular mass', 'temperature'],
        outputs=['mean molecular speed']
    )

    dsp.add_function(
        function=fun.calculate_molecular_diffusion_constant,
        inputs=['molecular mass', 'gas density', 'temperature', 'pressure'],
        outputs=['diffusion constant']
    )

    dsp.add_function(
        function=fun.calculate_particle_critical_diameter,
        inputs=['mean molecular speed', 'diffusion constant', 'accommodation '
                                                             'coefficient'],
        outputs=['particle critical diameter']
    )

    dsp.add_function(
        function=fun.calculate_particle_volume,
        inputs=['particle diameter'],
        outputs=['particle volume']
    )

    dsp.add_function(
        function=fun.calculate_particle_mass,
        inputs=['particle volume', 'gas density'],
        outputs=['particle mass']
    )

    dsp.add_function(
        function=fun.calculate_reduced_mass,
        inputs=['particle mass', 'molecular mass'],
        outputs=['reduced mass']
    )

    dsp.add_function(
       function=fun.calculate_center_mass_speed_1,
       inputs=['reduced mass', 'temperature'],
       outputs=['center of mass speed']
     )

    dsp.add_function(
        function=fun.calculate_speed_enhancement_factor,
        inputs=['particle mass', 'molecular mass'],
        outputs=['speed enhancement factor']
    )

    dsp.add_function(
        function=fun.calculate_center_mass_speed_2,
        inputs=['mean molecular speed', 'speed enhancement factor'],
        outputs=['center of mass speed']
    )

    dsp.add_function(
        function=fun.calculate_molecular_size_enhancement,
        inputs=['particle diameter', 'molecular diameter'],
        outputs=['molecular size enhancement']
    )

    dsp.add_function(
        function=fun.calculate_transition_regime_correction,
        inputs=['accommodation coefficient', 'knudsen number'],
        outputs=['transition correction']
    )

    dsp.add_function(
        function=fun.calculate_knudsen_number,
        inputs=['particle diameter', 'particle critical diameter'],
        outputs=['knudsen number']
    )

    dsp.add_function(
        function=fun.calculate_full_deposition_speed_1,
        inputs=['accommodation coefficient', 'molecular size enhancement',
               'transition correction', 'center of mass speed'],
        outputs=['full deposition speed']
    )

    dsp.add_function(
        function=fun.calculate_full_deposition_speed_2,
        inputs=['accommodation coefficient', 'molecular size enhancement',
                'speed enhancement factor',
                'mean molecular speed',
                'transition correction'],
        outputs=['full deposition speed']
    )

    dsp.add_function(
        function=fun.calculate_kinetic_deposition_speed,
        inputs=['mean molecular speed', 'accommodation coefficient'],
        outputs=['kinetic deposition speed']
    )

    dsp.add_function(
        function=fun.calculate_continuum_deposition_speed,
        inputs=['particle diameter', 'diffusion constant'],
        outputs=['continuum deposition speed']
    )

    dsp.add_function(
        function=fun.calculate_deposition_speed_normalised_ratio,
        inputs=['full deposition speed', 'kinetic deposition speed'],
        outputs=['deposition speed normalised ratio']
    )

    dsp.add_function(
        function=fun.calculate_diameter_growth_rate,
        inputs=['full deposition speed', 'gas density'],
        outputs=['diameter growth rate']
    )

    dsp.add_function(
        function=fun.calculate_kelvin_term,
        inputs=['particle diameter', 'kelvin diameter'],
        outputs=['kelvin term']
    )

    dsp.add_function(
        function=fun.calculate_evaporation_timescale,
        inputs=['particle diameter', 'kelvin term', 'full deposition '
                                                    'speed',
                'saturation ''concentration', 'gas density'],
        outputs=['evaporation timescale']
    )

    dsp.add_function(
        function=fun.calculate_particle_surface,
        inputs=['particle diameter'],
        outputs=['particle surface']
    )

    dsp.add_function(
        function=fun.calculate_collision_frequency_on_particles,
        inputs=['particle surface',
                'full deposition speed',
                'molecular mass', 'vapour concentration'],
        outputs=['collision frequency on particles']
    )

    # NEW FUNCTIONS
    dsp.add_function(
        function=fun.calculate_concentrations_Cs_i,
        inputs=['concentrations Cs_ip'],
        outputs=['concentrations Cs_i']
    )

    dsp.add_function(
        function=fun.calculate_concentrations_Cs_p,
        inputs=['concentrations Cs_ip'],
        outputs=['concentrations Cs_p']
    )

    dsp.add_function(
        function=fun.calculate_concentration_Cs_OA,
        inputs=['concentrations Cs_i'],
        outputs=['concentration Cs_OA']
    )

    dsp.add_function(
        function=fun.calculate_concentration_Cs_seed,
        inputs=['concentrations Cs_seed_p'],
        outputs=['concentration Cs_seed']
    )

    dsp.add_function(
        function=fun.calculate_number_concentration_Ns,
        inputs=['number concentrations Ns_p'],
        outputs=['number concentration Ns ']
    )

    dsp.add_function(
        function=fun.set_precursor_loss_Lprec,
        inputs=['production_time', 'precursor loss Lprec'],
        outputs=['precursor loss Lprec']
    )

    dsp.add_function(
        function=fun.calculate_net_vapour_productions_Pv_i,
        inputs=['fraction yields yi', 'concentrations Cs_p'],
        outputs=['net vapour productions Pv_i']
    )

    dsp.add_function(
        function=fun.calculate_activities_as_ip,
        inputs=['concentrations Cs_ip', 'precursor loss Lprec'],
        outputs=['activities as_ip']
    )

    dsp.add_function(
        function=fun.calculate_condensation_driving_forces_Fvs_ip,
        inputs=['activities as_ip', 'saturation concentrations Co_i',
                'concentrations Cv_i'],
        outputs=['condensation driving forces Fvs_ip']
    )

    dsp.add_function(
        function=fun.calculate_condensation_driving_forces_Fvs_ip,
        inputs=['activities as_ip', 'saturation concentrations Co_i',
                'concentrations Cv_i'],
        outputs=['condensation driving forces Fvs_ip']
    )

    dsp.add_function(
        function=fun.calculate_particle_volumes_vs_p,
        inputs=['number concentrations Ns_p',
                'seed density', 'organics density', 'concentrations Cs_p',
                'concentrations Cs_seed_p'],
        outputs=['particle volumes vs_p']
    )

    dsp.add_function(
        function=fun.calculate_particle_diameters_ds_p,
        inputs=['particle volumes vs_p'],
        outputs=['particle diameters ds_p']
    )

    dsp.add_function(
        function=fun.calculate_collision_frequencies_nus_ip,
        inputs=['full deposition speed', 'particle diameters ds_p',
                'number concentrations Ns_p'],
        outputs=['collision frequencies nus_ip']
    )

    dsp.add_function(
        function=fun.calculate_condensation_sinks_ks_ip,
        inputs=['collision frequencies nus_ip',
                'effective accommodation coefficient ai'],
        outputs=['condensation sinks ks_ip']
    )

    dsp.add_function(
        function=fun.calculate_condensation_fluxes_phi_vs_ip,
        inputs=['condensation sinks ks_ip',
                'condensation driving forces Fvs_ip'],
        outputs=['condensation fluxes phi_vs_ip']
    )

    dsp.add_function(
        function=fun.calculate_vapours_change_rates_dCv_i,
        inputs=['condensation fluxes phi_vs_ip',
                'net vapour productions Pv_i'],
        outputs=['vapours change rates dCv_i']
    )

    dsp.add_function(
        function=fun.calculate_suspended_particles_change_rates_dCs_ip,
        inputs=['condensation fluxes phi_vs_ip',
                'net particle productions Ps_ip'],
        outputs=['suspended particles change rates dCs_ip']
    )

    return dsp

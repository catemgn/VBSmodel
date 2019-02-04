import schedula as sh
import functions as fun


def model():
    """
    Returns the VBS Dynamics equations fof concentrations Cv_i and Cs_ip.
    WARNING! make sure that the order of the inputs is exactly the same of the
    order of inputs in the corresponding function.

    :return:
        workflow map of the model.
    :rtype: schedula.Dispatcher
    """

    # DEPOSITION SPEED MODULE

    dsp = sh.Dispatcher(
        name='calculate deposition speeds dynamics',
        description='Calculates the deposition speeds dynamics of the VBS.'
    )

    dsp.add_function(
        function=fun.calculate_mean_molecular_speed,
        inputs=['vapour masses m_i', 'temperature'],
        outputs=['mean molecular speed']
    )

    dsp.add_function(
        function=fun.calculate_molecular_diffusion_constant,
        inputs=['vapour masses m_i', 'organics density rho_org', 'temperature',
                'pressure'],
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
        inputs=['particle diameters ds_p'],
        outputs=['particle volume']
    )

    dsp.add_function(
        function=fun.calculate_particle_mass,
        inputs=['particle volume', 'organics density rho_org'],
        outputs=['particle mass']
    )

    dsp.add_function(
        function=fun.calculate_reduced_mass,
        inputs=['particle mass', 'vapour masses m_i'],
        outputs=['reduced mass']
    )

    dsp.add_function(
        function=fun.calculate_center_mass_speed,
        inputs=['reduced mass', 'temperature'],
        outputs=['center of mass speed']
    )

    dsp.add_function(
        function=fun.calculate_molecular_size_enhancement,
        inputs=['particle diameters ds_p', 'vapour effective diameters d_i'],
        outputs=['molecular size enhancement']
    )

    dsp.add_function(
        function=fun.calculate_transition_regime_correction,
        inputs=['accommodation coefficient', 'knudsen number'],
        outputs=['transition correction']
    )

    dsp.add_function(
        function=fun.calculate_knudsen_number,
        inputs=['particle diameters ds_p', 'particle critical diameter'],
        outputs=['knudsen number']
    )

    dsp.add_function(
        function=fun.calculate_full_deposition_speed,
        inputs=['accommodation coefficient', 'molecular size enhancement',
                'transition correction', 'center of mass speed'],
        outputs=['full deposition speeds']
    )

    dsp.add_function(
        function=fun.calculate_kinetic_deposition_speed,
        inputs=['mean molecular speed', 'accommodation coefficient'],
        outputs=['kinetic deposition speed']
    )

    dsp.add_function(
        function=fun.calculate_continuum_deposition_speed,
        inputs=['particle diameters ds_p', 'diffusion constant'],
        outputs=['continuum deposition speed']
    )

    dsp.add_function(
        function=fun.calculate_deposition_speed_normalised_ratio,
        inputs=['full deposition speeds', 'kinetic deposition speed'],
        outputs=['deposition speed normalised ratio']
    )

    dsp.add_function(
        function=fun.calculate_diameter_growth_rate,
        inputs=['full deposition speeds', 'organics density rho_org'],
        outputs=['diameter growth rate']
    )

    dsp.add_function(
        function=fun.calculate_evaporation_timescale,
        inputs=['particle diameters ds_p', 'kelvin terms Kp', 'full deposition '
                                                          'speeds',
                'saturation concentrations Co_i', 'organics density rho_org'],
        outputs=['evaporation timescale']
    )

    dsp.add_function(
        function=fun.calculate_particle_surface,
        inputs=['particle diameters ds_p'],
        outputs=['particle surface']
    )

    dsp.add_function(
        function=fun.calculate_collision_frequency_on_particles,
        inputs=['particle surface',
                'full deposition speeds',
                'vapour masses m_i', 'concentrations Cv_i'],
        outputs=['collision frequency on particles']
    )

    # DYNAMICS EQUATIONS MODULE

    dsp2 = sh.Dispatcher(
        name='VBS dynamics equations',
        description='Provides the VBS dynamics equations.'
    )

    dsp.add_data(data_id='inputs')
    dsp2.add_data(data_id='number of species ni')
    dsp2.add_data(data_id='number of populations np')
    dsp2.add_data(data_id='saturation concentrations Co_i')
    dsp2.add_data(data_id='concentrations Cs_ip')
    dsp2.add_data(data_id='concentrations Cv_i')
    dsp2.add_data(data_id='concentrations Cs_ip')
    dsp2.add_data(data_id='net particle productions Ps_ip')
    dsp2.add_data(data_id='fraction yields yi')
    dsp2.add_data(data_id='integration time')
    dsp2.add_data(data_id='precursor loss Lprec')
    dsp2.add_data(data_id='number concentrations Ns_p')
    dsp2.add_data(data_id='organics density rho_org')
    dsp2.add_data(data_id='seed density rho_seed_p')
    dsp2.add_data(data_id='seed diameters d_seed_p')
    dsp2.add_data(data_id='number concentrations Ns_p')
    dsp2.add_data(data_id='vapour effective diameters d_i')
    dsp2.add_data(data_id='vapour masses m_i')
    dsp2.add_data(data_id='temperature', default_value=298)  # temperature [K].
    dsp2.add_data(data_id='pressure', default_value=101325)  # pressure [Pa].
    dsp2.add_data(data_id='kelvin diameter', default_value=4.5) # kelvin # diameter for organics [nm].
    dsp2.add_data(data_id='accommodation coefficient', defaults_value=1)
    dsp2.add_data(data_id='effective accommodation coefficient',
                  defaults_value=1)

    dsp2.add_dispatcher(
        dsp,
        dsp_id='calculate deposition speeds dynamics',
        inputs=['particle diameters ds_p', 'organics density rho_org',
                'vapour effective diameters d_i', 'vapour masses m_i',
                'temperature', 'saturation concentrations Co_i', 'pressure',
                'kelvin terms Kp', 'accommodation coefficient', 'concentrations Cv_i'],
        outputs=['full deposition speeds', 'collision frequency on '
                                           'particles', 'evaporation timescale']
    )

    dsp2.add_function(
        function=fun.calculate_concentrations_Cs_i,
        inputs=['concentrations Cs_ip'],
        outputs=['concentrations Cs_i']
    )

    dsp2.add_function(
        function=fun.calculate_concentrations_Cs_p,
        inputs=['concentrations Cs_ip'],
        outputs=['concentrations Cs_p']
    )

    dsp2.add_function(
        function=fun.calculate_concentration_Cs_OA,
        inputs=['concentrations Cs_i'],
        outputs=['concentration Cs_OA']
    )

    dsp2.add_function(
        function=fun.calculate_concentrations_Cs_seed_p,
        inputs=['seed density rho_seed_p', 'seed diameters d_seed_p',
                 'number concentrations Ns_p'],
        outputs=['concentrations Cs_seed_p']
    )

    dsp2.add_function(
        function=fun.calculate_concentration_Cs_seed,
        inputs=['concentrations Cs_seed_p'],
        outputs=['concentration Cs_seed']
    )

    dsp2.add_function(
        function=fun.calculate_number_concentration_Ns,
        inputs=['number concentrations Ns_p'],
        outputs=['number concentration Ns ']
    )

    dsp2.add_function(
        function=fun.set_precursor_loss_Lprec,
        inputs=['integration time', 'precursor loss Lprec'],
        outputs=['precursor loss Lprec']
    )

    dsp2.add_function(
        function=fun.calculate_net_vapour_productions_Pv_i,
        inputs=['fraction yields yi', 'precursor loss Lprec'],
        outputs=['net vapour productions Pv_i']
    )

    dsp2.add_function(
        function=fun.calculate_particle_volumes_vs_p,
        inputs=['number concentrations Ns_p',
                'seed density rho_seed_p', 'organics density rho_org',
                'concentrations Cs_p',
                'concentrations Cs_seed_p'],
        outputs=['particle volumes vs_p']
    )

    dsp2.add_function(
        function=fun.calculate_particle_diameters_ds_p,
        inputs=['particle volumes vs_p'],
        outputs=['particle diameters ds_p']
    )

    dsp2.add_function(
        function=fun.calculate_kelvin_terms_Kp,
        inputs=['particle diameters ds_p', 'kelvin diameter'],
        outputs=['kelvin terms Kp']
    )

    dsp2.add_function(
        function=fun.calculate_collision_frequencies_nus_ip,
        inputs=['full deposition speeds', 'particle diameters ds_p',
                'number concentrations Ns_p'],
        outputs=['collision frequencies nus_ip']
    )

    dsp2.add_function(
        function=fun.calculate_condensation_sinks_ks_ip,
        inputs=['collision frequencies nus_ip',
                'effective accommodation coefficient'],
        outputs=['condensation sinks ks_ip']
    )

    dsp2.add_function(
        function=fun.calculate_activities_as_ip,
        inputs=['concentrations Cs_ip', 'concentrations Cs_p'],
        outputs=['activities as_ip']
    )

    dsp2.add_function(
        function=fun.calculate_condensation_driving_forces_Fvs_ip,
        inputs=['activities as_ip', 'saturation concentrations Co_i',
                'concentrations Cv_i', 'kelvin terms Kp'],
        outputs=['condensation driving forces Fvs_ip']
    )

    dsp2.add_function(
        function=fun.calculate_condensation_fluxes_phi_vs_ip,
        inputs=['condensation sinks ks_ip',
                'condensation driving forces Fvs_ip'],
        outputs=['condensation fluxes phi_vs_ip']
    )

    dsp2.add_function(
        function=fun.calculate_vapours_change_rates_dCv_i,
        inputs=['condensation fluxes phi_vs_ip',
                'net vapour productions Pv_i'],
        outputs=['vapours change rates dCv_i']
    )

    dsp2.add_function(
        function=fun.calculate_suspended_particles_change_rates_dCs_ip,
        inputs=['condensation fluxes phi_vs_ip',
                'net particle productions Ps_ip'],
        outputs=['suspended particles change rates dCs_ip']
    )

    return dsp2


# def define_differential_equations(inputs):
#     """"
#     Return the VBS Dynamics equations in a form to be passed to the ode solver.
#     """
#     dsp = model.copy()
#     for k, v in dsp(defaults).items():
#         dsp.set_default_value(k, v)
#     dsp.add_function(
#          ‘split_inputs’,
#           lambda i, j, x: (x[:i], x[i:].reshape(i,j)),
#           [‘ni’, ‘np’, ‘inputs’],  # N.b. ‘ni’ and ‘np’ have to be defined into defaults.
#           [‘cvi’, ‘csip’]
#     )
#     dsp.add_function(
#          ‘merge_outputs’,
#           lambda x, y: np.append(x, y.ravel()),
#           [‘dcvi’, ‘dcsip’],
#           [‘outputs’],
#     )
#     return sh.DispatchPipe(dsp, [‘inputs’], [‘outputs’])



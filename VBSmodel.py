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
        function=fun.calculate_collision_frequency,
        inputs=['particle surface',
                'full deposition speed',
                'molecular mass', 'vapour concentration'],
        outputs=['collision frequency']
    )

    return dsp

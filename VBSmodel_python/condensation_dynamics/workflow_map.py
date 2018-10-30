import schedula as sh
from VBSmodel_python.condensation_dynamics import functions as fun


def model():
    """
    Returns a workflow map of the model (dispatcher).

    :return:
        workflow map of the model.
    :rtype: schedula.Dispatcher
    """

    dsp = sh.Dispatcher(
        name='deposition-speed',
        description='Model that calculates the deposition speed of a vapour '
                    'oer a particle.'
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
        inputs=['mean molecular speed', 'diffusion constant'],
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
        inputs=['particle diameter', 'particle critical diameter',
                'accommodation coefficient'],
        outputs=['transition correction']
    )

    dsp.add_function(
        function=fun.calculate_knudsen_number,
        inputs=['particle diameter', 'particle critical diameter'],
        outputs=['knudsen_number']
    )

    dsp.add_function(
        function=fun.calculate_full_deposition_speed,
        inputs=['accommodation coefficient', 'molecular size enhancement',
                'transition correction', 'center of mass speed'],
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
                'speed', 'saturation ''concentration', 'gas density'],
        outputs=['evaporation timescale']
    )
    return dsp

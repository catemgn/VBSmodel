dfl = {
    'values': {
        # generic parameters/default values used in the code.

        'temperature': 300,  # reference ambient temperature [K].
        'temperature2': 298,  # other reference temperature used [K].
        'pressure': 101325,  # 1 atm pressure [Pa].
        'gas_density': 1.4e3,  # typical organics density [Kg/m3].
        'gas_molar_weight': 250,  # typical organics gas molar weight [g/mol].
        'air_mass': 29,  # mass of air [g/mol].
        'air_density': 874,  # density of liquid air [kg/m3].
        'molecular_diameter': 0.9,  # typical organic gas diameter [nm].
        'C*': 1  # saturation concentration [ug/m3].
        },


    'functions': {
        # specific parameters used in some functions, listed per function name.

        'molecular_diffusion_constant': {
            'm': 1.7,  # temperature_scaling_exponent due to the temperature
            # dependence of the collision integral [-].
            'delta': 2/3,  # _scaling_factor that  relates the bulk effective
            # diameter and the collisional diameter [-].
        },
        'transition_regime_correction': {
            'dp_critic': 110,  # particle critical diameter [nm].
            'alpha': 1  # accommodation parameter [-].
        }
    },

    'conversions': {
        # useful conversion factors

        'kg2g': 1e3,  # kilograms to grams.
        'kg2ug': 1e9,  # kilograms to micrograms.
        'nm2m': 1e-9,  # nanometer to meters.
        'kg2kg/mol': 1.660538921e-27,  # kg to kg/mol.
        'sec2hr': 1/3.6e3  # sec to hours.

    }
}

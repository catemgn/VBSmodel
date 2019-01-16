import numpy as np
from VBSmodel_python.condensation_dynamics.workflow_map import model

# settings variables

particle_diameter = np.arange(5e-1, 5e4, 1e2)
molecular_mass = 250
molecular_density = 1400
molecular_diameter = 0.8
temperature = 298
pressure = 101326
accom_coeff = 1

inputs = {'molecular mass': molecular_mass, 'gas density': molecular_density,
          'temperature': temperature, 'pressure':pressure}

sol_dsp = model().dispatch(inputs=inputs, outputs={'particle critical '
                                                   'diameter'})
print(sol_dsp)
sol_dsp.plot()

#'particle diameter': particle_diameter, 'molecular_diameter': molecular_diameter,

# inputs={'temperature': temperature, 'accommodation coefficient': accom_coeff,
#              'molecular mass': molecular_mass}

# inputs2={'temperature': temperature, 'accommodation coefficient': accom_coeff,
#              'molecular mass': molecular_mass}
# sol_dsp = model().dispatch(inputs=inputs2, outputs=['kinetic deposition speed'])
# print(sol_dsp)
# sol_dsp.plot()
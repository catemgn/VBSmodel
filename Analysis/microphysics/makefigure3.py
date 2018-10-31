import numpy as np
from VBSmodel_python.condensation_dynamics.workflow_map import model

# settings variables

particle_diameter = np.arange(1, 100, 1)  #np.arange(5e-1, 5e4, 1e2)
molecular_mass = 250
molecular_density = 1400
molecular_diameter = 0.9
temperature = 298
pressure = 101326
accom_coeff = 1
kelvin_d = 4.5
c = 1
cv=1

inputs = {'particle diameter': particle_diameter, 'molecular mass':
     molecular_mass, 'gas density': molecular_density,
          'molecular diameter': molecular_diameter, 'accommodation '
                                                    'coefficient': accom_coeff,
          'temperature': temperature, 'pressure': pressure, 'kelvin '
                                                            'diameter':
              kelvin_d, 'saturation concentration': c, 'vapour '
                                                       'concentration':cv}


sol_dsp = model().dispatch(inputs=inputs, outputs={'collision frequency'})
print(sol_dsp)
sol_dsp.plot()


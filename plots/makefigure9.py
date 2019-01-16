""""
Plot figure 9 paper.
Direct use of the VBS model (no process).
"""

import numpy as np
import matplotlib.pyplot as plt
from VBSmodel import model

# settings variables
particle_diameter = np.logspace(-0.5, 4.5, num=501, base=10, endpoint=True)

inputs = {'particle diameter': particle_diameter,
          'molecular mass': 250, 'gas density': 1400,
          'molecular diameter': 0.8,
          'accommodation coefficient': 1,
          'temperature': 298,
          'pressure': 101326}

# PLOT FOR alpha=1
sol_dsp = model().dispatch(inputs=inputs)  # Calculating model solution for
# the given inputs.
plt.loglog(particle_diameter, sol_dsp['full deposition speed'], '-', c='g',
           linewidth=3, label='Organics alpha=1')
# Plot asymptotes
selector = (particle_diameter < sol_dsp['particle critical diameter'])
dp = particle_diameter[selector]
v_kinetic = np.full(len(dp), sol_dsp['kinetic deposition speed'])
plt.loglog(dp, v_kinetic, '--', c='g')
dp = particle_diameter[~selector]
v_continuum = sol_dsp['continuum deposition speed'][~selector]
plt.loglog(dp, v_continuum, '--', c='g')

# PLOT FOR alpha=0.1
inputs['accommodation coefficient'] = 0.1
sol_dsp = model().dispatch(inputs=inputs)
plt.loglog(particle_diameter, sol_dsp['full deposition speed'], '-.',
           c='g', linewidth=2, label='Organics alpha=0.1')
# Plot asymptotes
selector = (particle_diameter < sol_dsp['particle critical diameter'])
dp = particle_diameter[selector]
v_kinetic = np.full(len(dp), sol_dsp['kinetic deposition speed'])
plt.loglog(dp, v_kinetic, '-.', c='g')
dp = particle_diameter[~selector]
v_continuum = sol_dsp['continuum deposition speed'][~selector]
plt.loglog(dp, v_continuum, '-.', c='g')


# PLOT OFR NITRIC ACID
accommodation_coefficient = 1
molecular_mass = 63
molecular_density = 1700
molecular_diameter = 0.4
inputs['accommodation coefficient'] = 1
inputs['molecular mass'] = 63
inputs['gas density'] = 1700
inputs['molecular diameter'] = 0.4
sol_dsp = model().dispatch(inputs=inputs)
plt.loglog(particle_diameter, sol_dsp['full deposition speed'], '-', c='b',
           linewidth=2, label='Nitric acid')
# Plot asymptotes
selector = (particle_diameter < sol_dsp['particle critical diameter'])
dp = particle_diameter[selector]
v_kinetic = np.full(len(dp), sol_dsp['kinetic deposition speed'])
plt.loglog(dp, v_kinetic, '--', c='b')
dp = particle_diameter[~selector]
v_continuum = sol_dsp['continuum deposition speed'][~selector]
plt.loglog(dp, v_continuum, '--', c='b')

# Setting figure parameters
plt.ylim((1, 1e3))
plt.xlim((.3, 3e4))
plt.xlabel('particle diameter d_p(nm)')
plt.ylabel('deposition speed (m/s)')
plt.text(.7, 550, 'Collision', fontsize=11)
plt.text(10, 120, 'Kinetic', fontsize=11)
plt.text(2000, 25, 'Continuum', fontsize=11)
plt.legend(loc='best')
plt.show()


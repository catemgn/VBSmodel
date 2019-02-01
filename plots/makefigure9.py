""""
Plot figure 9 paper.
Direct use of the VBSDynamics model (no process).
"""

import numpy as np
import matplotlib.pyplot as plt
from VBSmodel import model

# settings variables
particle_diameter = np.logspace(-0.5, 4.5, num=501, base=10, endpoint=True)

inputs = {'particle diameters ds_p': particle_diameter,
          'vapour masses m_i': np.array([250]),
          'organics density rho_org': 1400,
          'vapour effective diameters d_i': np.array([0.8]),
          'accommodation coefficient': 1,
          'temperature': 298,
          'pressure': 101326,
          'kelvin diameter': 4.5}

dep_sp = 'full deposition speeds'
crit_d = 'particle critical diameter'
kin_sp = 'kinetic deposition speed'
con_sp = 'continuum deposition speed'


# PLOT FOR alpha=1
sol_dsp = model().dispatch(inputs=inputs)  # get the model solution
# given inputs.
plt.loglog(particle_diameter, sol_dsp[dep_sp].flatten(), '-', c='g',
           linewidth=3,
           label='Organics alpha=1')

# Plot asymptotes
selector = (particle_diameter < sol_dsp.get_node('calculate deposition '
                                                 'speeds', crit_d)[0])
dp = particle_diameter[selector]
v_kinetic = np.full(len(dp), sol_dsp.get_node(
            'calculate deposition speeds', kin_sp)[0])
plt.loglog(dp, v_kinetic, '--', c='g')
dp = particle_diameter[~selector]
v_continuum = sol_dsp.get_node(
    'calculate deposition speeds', con_sp)[0].flatten()[~selector]
plt.loglog(dp, v_continuum, '--', c='g')


# PLOT FOR alpha=0.1
inputs['accommodation coefficient'] = 0.1
sol_dsp = model().dispatch(inputs=inputs)
plt.loglog(particle_diameter, sol_dsp[dep_sp].flatten(), '-.',
           c='g', linewidth=2, label='Organics alpha=0.1')

# Plot asymptotes
selector = (particle_diameter < sol_dsp.get_node('calculate deposition '
                                                 'speeds', crit_d)[0])
dp = particle_diameter[selector]
v_kinetic = np.full(len(dp),
                    sol_dsp.get_node('calculate deposition speeds', kin_sp)[0])
plt.loglog(dp, v_kinetic, '-.', c='g')
dp = particle_diameter[~selector]
v_continuum = sol_dsp.get_node('calculate deposition speeds', con_sp)[
    0].flatten()[~selector]
plt.loglog(dp, v_continuum, '-.', c='g')


# PLOT OFR NITRIC ACID
inputs['accommodation coefficient'] = 1
inputs['vapour masses m_i'] = np.array([63])
inputs['organics density rho_org'] = 1700
inputs['vapour effective diameters d_i'] = np.array([0.4])
sol_dsp = model().dispatch(inputs=inputs)
plt.loglog(particle_diameter, sol_dsp[dep_sp].flatten(), '-', c='b',
           linewidth=2, label='Nitric acid')

# Plot asymptotes
selector = (particle_diameter < sol_dsp.get_node('calculate deposition '
                                                 'speeds', crit_d)[0])
dp = particle_diameter[selector]
v_kinetic = np.full(len(dp),
                    sol_dsp.get_node('calculate deposition speeds', kin_sp)[0])
plt.loglog(dp, v_kinetic, '--', c='b')
dp = particle_diameter[~selector]
v_continuum = sol_dsp.get_node('calculate deposition speeds', con_sp)[
    0].flatten()[~selector]
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

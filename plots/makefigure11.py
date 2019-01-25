""""
Plot figure 10 paper.
Direct use of the VBS model (no process).
"""

import numpy as np
import matplotlib.pyplot as plt
from VBSmodel import model


def plot_collision_frequency(inputs, cv):
    """
    Plot the frequency with at which molecules collide with particles as a
    function of particle diameter for given vapor concentrations (Cv) [ug/m3].
    """
    inputs['vapour concentration'] = cv
    sol_dsp = model().dispatch(inputs=inputs, outputs=['collision frequency'])
    mu = sol_dsp['collision frequency']

    return plt.loglog(inputs['particle diameter'], mu, '--', c='k')


inputs = {'particle diameter': np.logspace(-0.5, 4.5, num=501, base=10,
                                           endpoint=True),
          'molecular mass': 250,
          'gas density': 1400,
          'molecular diameter': 0.9,
          'accommodation coefficient': 1,
          'temperature': 298,
          'pressure': 101326
          }

p10 = plot_collision_frequency(inputs, 1)
plt.setp(p10, ls='-', linewidth=3, color=(0, 0.4, 0), label='Cv = 1 ug/m3')
p1m2 = plot_collision_frequency(inputs, 1e-2)
plt.setp(p1m2, ls='-.', linewidth=2, color=(0, 0.2, 0), label='Cv = 0.01 '
                                                              'ug/m3')
plot_collision_frequency(inputs, 0.1)
plot_collision_frequency(inputs, 1e-3)

# Setting figure parameters
plt.xlim((0.3, 3e4))
plt.ylim((3e-3, 3e5))
plt.xlabel('particle diameter d_p(nm)')
plt.ylabel('Collision frequency for particles  (1/s)')
plt.legend(loc='best', fontsize=11)
plt.grid(True, which="both", ls="--", alpha=0.5)

plt.show()

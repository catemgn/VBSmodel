""""
Plot figure 10 paper.
Direct use of the VBSDynamics model (no process).
"""

import numpy as np
import matplotlib.pyplot as plt
from VBSmodel import model


def plot_evaporation_timescale(inputs, co):
    """
    plots the timescale for evaporation of a constituent with a given
    saturation concentration co [ug/m3] from particles of various sizes.
    """

    inputs['saturation concentrations Co_i'] = co

    sol_dsp = model().dispatch(inputs=inputs)
    tau = sol_dsp.get_node('calculate deposition speeds', 'evaporation '
                                                          'timescale')[0]

    return plt.loglog(inputs['particle diameters ds_p'], tau, '--', c='k')


inputs = {'particle diameters ds_p': np.logspace(-0.5, 4.5, num=501, base=10,
                                           endpoint=True),
          'vapour masses m_i': 350,
          'organics density rho_org': 1400,
          'vapour effective diameters d_i': 0.8,
          'accommodation coefficient': 1,
          'temperature': 300,
          'pressure': 101325,
          'kelvin diameter': 4.5
          }

p1m2 = plot_evaporation_timescale(inputs, 1e-2)
plt.setp(p1m2, ls='-.', linewidth=2, color=(0, 0.2, 0), label='C^o = 0.01 '
                                                              'ug/m3')
p10 = plot_evaporation_timescale(inputs, 1)
plt.setp(p10, ls='-', linewidth=3, color=(0, 0.4, 0), label='C^o =1 '
                                                            'ug/m3')
p12 = plot_evaporation_timescale(inputs, 1e2)
plt.setp(p12, ls='-.', linewidth=2, color=(0, 0.6, 0), label='C^o = 100 '
                                                             'ug/m3')

p14 = plot_evaporation_timescale(inputs, 1e4)
plt.setp(p14, ls='-.', linewidth=2, color=(0, 0.8, 0), label='C^o = 10^4 '
                                                             'ug/m3')

plot_evaporation_timescale(inputs, 1e-3)
plot_evaporation_timescale(inputs, 0.1)
plot_evaporation_timescale(inputs, 10)
plot_evaporation_timescale(inputs, 1e3)

# Setting figure parameters
plt.xlim((0.3, 3e4))
plt.ylim((1e-3, 1e3))
plt.xlabel('particle diameter d_p(nm)')
plt.ylabel('evaporation timescale (hr)')
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend(loc='best', fontsize=9)

plt.show()


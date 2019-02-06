import numpy as np

"""
This dictionary provides the input setup for calculations.
"""

VBSsetup = {

    # general properties:
  'parameters':{

    'T': 300,  # temperature [K}.
    'P': 101325,  # pressure [Pa].
    'nPop': 1,  # number of the seed particle populations [-].
    'nOrg': 4,  # number of species [-].
    'Co_i': np.array([0.01, 1, 100, 10000]),  # saturation concentrations set (Volatility set) distribution [ug/m3].
    'y': np.array([0.25, 0.25, 0.25, 0.25]),  # yields of products. Add to 1 [-].
    'Lprec': 1,  # precursor loss rate [ug/(m3*min)]. # TODO min or h?

    # seed and organic properties

    'diamSeed': np.array([100]),  # diameter of the  seed particle in the
    # populations [nm].
    'rhoSeed': np.array([1770]),  # densities of the seed particle of
    # population [kg/m3].
    'Ns_p': np.array([16000]),  # number concentration of suspended
    # particles for each populations [1/cm3].
    'mOrg': np.ones(4)*350,  # molar weight of organics [g/mol].
    'diamOrg': np.ones(4)*0.8,  # effective vapour diameter for each
    # species# TODO units?.
    'alphaOrg': 1,  # mass accommodation coefficient for organics [-].
    'rhoOrg': 1400,  # density of organics [kg/m3].
    'diamK10org': 4.5,  # kelvin diameter for organics [nm].
    'net particle productions Ps_ip': np.full((4, 1), 0),
    'effectiveai': 1
},

    'initial conditions': {

    'time': (-100, 0),
    # initial conditions
    'Cv_i_init': np.zeros(4),  # initial vapours
    # concentrations for each
    # species i [ug/m3].
    'Cs_ip_init': np.full((4, 1), 0)  # initial suspended
    # particle
    # concentrations for each species i (row) and population p (column) [ug/m3].
}
}


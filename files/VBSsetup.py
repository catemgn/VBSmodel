import numpy as np

"""
This dictionary provides a default setup for calculations.
Note that most arrays need to be columns for differential equations solvers.
"""

"""
VBSsetup = {

    # general properties:

    'T': 300,  # temperature [K}.
    'P': 101325,  # pressure [Pa].
    'nOrg': 4,  # number of organic species [-].
    'd': np.array([0.01, 1, 100, 10000]),  # volatility set distribution [
    # ug/m3].
    'y': np.array([[0.25, 0.25, 0.25, 0.25]]),  # yields of products. Add to 1
    #  [-].
    'LPrec': 1,  # reaction rate [ug/(m3*min)].

    # properties of the seed particle population:

    'nPop': 1,  # number of the seed particle population [-].
    'diamSeed': 100,  # diameter of the  seed particle in the population [nm].
    'rhoSeed': 1.77e3,  # density of the seed particle of population [kg/m3].
    'nSeed': 16000,  # number of seed particles in the population [1/cm3].

    # properties of the organics:

    'mOrg': 350,  # molar weight of organics [g/mol].
    'diamOrg': 0.8,  # effective vapour diameter TODO units?.
    'alphaOrg': 1,  # mass accommodation coefficient for organics [-].
    'rhoOrg': 1.4e3,  # density of organics [kg/m3].
    'diamK10org': 4.5,  # kelvin diameter for organics [nm].
"""

{
    #specific figure inputs #Figure3

    'particle diameter': np.logspace(-0.5, 4.5, num=501, base=10,
                                      endpoint=True),
    'molecular mass': 250,
    'molecular density': 1400,
    'molecular diameter': 0.8,
    'temperature': 298,
    'pressure': 101326,
    'accommodation coefficient': 1

}

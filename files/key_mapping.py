"""
This dictionary contains the map between the keys values of the input
dictionary and the keys values of the VBS model.
(Mapping order:  input_name : model_name)
"""

key_mapping = {
    'T': 'temperature',
    'P': 'pressure',
    'Co_i': 'saturation concentrations Co_i',
    'y': 'fraction yields yi',
    'Lprec': 'precursor loss Lprec',

    'rhoSeed': 'seed density rho_seed_p',
    'diamSeed': 'seed diameters d_seed_p',
    'Ns_p': 'number concentrations Ns_p',
    'mOrg': 'vapour masses m_i',
    'diamOrg': 'vapour effective diameters d_i',
    'rhoOrg': 'organics density rho_org',
    'alphaOrg': 'accommodation coefficient',
    'diamK10org': 'kelvin diameter',
    
    'Cs_ip_init': 'concentrations Cs_ip',
    'Cv_i_init': 'concentrations Cv_i',




    # figure 3 specific inputs

    # 'particle diameter': 'particle diameters ds_p',
    # 'molecular mass': 'vapour masses m_i',
    # 'molecular density': 'organics density rho_org',
    # 'molecular diameter': 'vapour effective diameters d_i',
    # 'accommodation coefficient': 'accommodation coefficient',



}
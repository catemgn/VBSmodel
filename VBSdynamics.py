import schedula as sh
import numpy as np
from VBSmodel import model
from scipy.integrate import odeint


def define_differential_equations(parameters):
    """"
    Return the VBS Dynamics equations in a form to be passed to the ode solver.
    """
    dsp = model().copy()

    # set the input parameters to be used in the model.
    for k, v in dsp(parameters).items():
        dsp.set_default_value(k, v)

    # split the inputs.
    dsp.add_function('split_initial_conditions',
                     lambda i, j, x: (x[:i], x[i:].reshape(i, j)),
                     ['number of species ni', 'number of populations np','inputs'],
                     ['concentrations Cv_i', 'concentrations Cs_ip']
                     )
    # pack outputs together.
    dsp.add_function('merge_outputs',
                     lambda x, y: np.append(x, y.ravel()),
                     ['vapours change rates dCv_i',
                      'suspended particles change rates dCs_ip'], ['outputs'])

    return sh.DispatchPipe(dsp, ['inputs'], ['outputs'])


def ode_solve_dynamics(concentrations, time, dyn_function):
    """
    Solve the equations of the VBS dynamics for the vapours and
    suspended particle concentrations.

    :param concentrations: vector of all concentrations
    :type: concentrations: np.array (1D).

    :param time: time period for integration.

    :param dyn_function: np.array (1D).

    :return: VBS dynamics equation for concentrations.
    :rtype: np.array
    """
    return odeint(dyn_function, concentrations, time)


def solve_dyn():

    solver = sh.Dispatcher(
        name='Solve VBS Dynamics',
        description='Solve the VBS Dynamics equations fof concentrations Cv_i and Cs_ip.'
    )

    solver.add_function('pack inputs',
                     lambda x, y: np.concatenate((x, y), axis=None),
                        ['concentrations Cv_i', 'concentrations Cs_ip'], ['inputs'])

    solver.add_data(data_id='parameters', default_value={})

    solver.add_function(
        function_id='define_differential_equations',
        function=define_differential_equations,
        inputs=['parameters'],
        outputs=['dynamics equations']
    )

    solver.add_function(
        function_id='solve VBS dynamics',
        function=ode_solve_dynamics,
        inputs=['inputs', 'integration time', 'dynamics equations'],
        outputs=['outputs']
    )

    return solver


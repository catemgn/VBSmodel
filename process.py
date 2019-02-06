"""
Defines the file processing the VBS Dynamics model.
"""

import schedula as sh
from VBSdynamics import solve_dyn
import pandas as pd

process = sh.BlueDispatcher(name='Processing Model')


# BluePrint create a 'canvas' for the Dispatcher. It is not compiled yet.


@sh.add_function(process, outputs=['raw inputs'])
def read_data(input_fpath):
    """
    Reads the input files

    :param input_fpath:
        Input file path.
    :type input_fpath: str

    :return:
        Raw Data.
    :rtype: dict
    """

    from files.VBSinputs import VBSsetup

    return VBSsetup


process.add_data(
    data_id='key_mapping', default_value={},
    description='Dictionary renaming keys.'
)


# TODO Find a more elegant and compact way to do this function.
@sh.add_function(process, inputs=['key_mapping', 'raw inputs'], outputs=[
    'inputs'])
def parse_data(nested_keys, nested_inputs):
    par = sh.map_dict(nested_keys['parameters'], nested_inputs['parameters'])
    ic = sh.map_dict(nested_keys['initial conditions'],
                     nested_inputs['initial conditions'])

    return sh.map_list(('parameters', 'initial conditions'), par, ic)


process.add_function(
    function_id='compute VBS Dynamics',
    function=sh.SubDispatch(solve_dyn()),
    inputs=['inputs'],
    outputs=['outputs'],
    description='Executes the computational model.'
)


@sh.add_function(process)
def save_outputs(outputs, output_fpath):
    """
    Save model outputs in a json dictionary.

    :param outputs:
        Model outputs.
    :type outputs: dict

    :param output_fpath:
        Output file path.
    :type output_fpath: str
    """

    df = pd.DataFrame(outputs)
    with pd.ExcelWriter(output_fpath) as writer:
        df.to_excel(writer)


if __name__ == '__main__':
    process.register().plot()  # register compiles the dispatcher.

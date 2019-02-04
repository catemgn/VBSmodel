"""
Defines the file processing the VBS Dynamics model.
"""

import schedula as sh
from VBSmodel import model
import pandas as pd

process = sh.BlueDispatcher(name='Processing Model')


# BluePrint create a 'canvas' for the Dispatcher. It is not compiled yet.


@sh.add_function(process, outputs=['raw parameters', 'raw initial conditions'])
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

    from files.VBSinputs import VBSsetup, dynamics_IC

    return VBSsetup, dynamics_IC


process.add_data(
    data_id='key_mapping', default_value={},
    description='Dictionary renaming keys.'
)

process.add_function('parse data',
                     lambda keys, dict1, dict2: (
                         sh.map_dict(keys, dict1), sh.map_dict(keys, dict2)),
                     ['key_mapping', 'raw parameters',
                      'raw initial conditions'],
                     ['parameters', 'initial conditions']
                     )

process.add_data(data_id='parameters')

process.add_function(
    function_id='compute VBS Dynamics',
    function=sh.SubDispatch(model()),
    inputs=['parameters', 'initial conditions'],
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

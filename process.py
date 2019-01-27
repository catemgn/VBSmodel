
"""
Defines the file processing the VBS model.
"""

import schedula as sh
import json
from VBSmodel import model


process = sh.BlueDispatcher(name='Processing Model')
# BluePrint create a 'canvas' for the Dispatcher. It is not compiled yet.


@sh.add_function(process, outputs=['raw_data'])
def read_data(input_fpath):
    """
    Reads the json file. Input data are a json dictionary.

    :param input_fpath:
        Input file path.
    :type input_fpath: str

    :return:
        Raw Data.
    :rtype: dict
    """
    with open(input_fpath, 'r') as writer:
        raw_data = json.load(writer)
    return raw_data


process.add_data(
    data_id='key_mapping', default_value={},
    description='Dictionary renaming keys.'
)

process.add_function(
    function_id='parse_data',
    function=sh.map_dict,
    inputs=['key_mapping', 'raw_data'],
    outputs=['data'],
    description='Rename the raw data names.'
)

process.add_function(
    function_id='execute_VBS',
    function=sh.SubDispatch(model()),
    inputs=['data'],
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
    with open(output_fpath, 'w') as writer:
        json.dump(outputs, writer)

if __name__ == '__main__':
    process.register().plot()   # register compiles the dispatcher.
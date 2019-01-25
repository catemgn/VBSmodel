"""
Defines the file processing the VBS model.
"""

import schedula as sh
import pandas as pd
from VBSmodel import model


process = sh.BlueDispatcher(name='Processing Model')
# BluePrint create a 'canvas' for the Dispatcher. It is not compiled yet.


@sh.add_function(process, outputs=['raw_data'])
def read_data(): #TODO change input read for different input files format.
    """
    Reads the excel file.

    :param input_fpath:
        Input file path.
    :type input_fpath: str

    :return:
        Raw Data.
    :rtype: dict
    """
    from files.VBSsetup import VBSsetup
    return VBSsetup


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
    Save model outputs in an Excel file.

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
    process.register().plot()   # register compiles the dispatcher.
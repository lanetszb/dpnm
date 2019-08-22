import sys
import os
import pandas as pd

from input.parameters_frame import ParametersFrame

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))


class Parameters(ParametersFrame):
    def __init__(s, config_file):
        super().__init__(config_file)

    def process_throats(s):
        s.throats = pd.read_csv(s.pore_throats, index_col=0)
        print(s.throats)

    def process_pores(s):
        s.pores = pd.read_csv(s.pores_data, index_col=0)
        print(s.pores)


if __name__ == '__main__':
    pore_throats = Parameters(config_file=sys.argv[1])
    pore_throats.process_throats()
    pore_throats.process_pores()

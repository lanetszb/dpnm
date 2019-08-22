import sys
import os
import configparser

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))


class ParametersFrame:
    def __init__(s, config_file):
        s.config = configparser.ConfigParser()
        s.config.read(config_file)

        s.pore_throats = str(s.config.get('PNData', 'pore_throats'))
        s.pores_data = str(s.config.get('PNData', 'pores_data'))

        s.throats = None
        s.pores = None

    def __str__(s):
        out_str = ' '
        if s.throats is not None:
            out_str = 'pore_throats ' + str(s.pore_throats)
        if s.pores is not None:
            out_str = 'pores_data ' + str(s.pores_data)
        return out_str


if __name__ == '__main__':
    parameters_frame = ParametersFrame(config_file=sys.argv[1])
    print(parameters_frame)

import sys
import os
import configparser
import numpy as np

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))


class Props:
    def __init__(s, config_file):
        s.__config = configparser.ConfigParser()
        s.__config.read(config_file)
        get = s.__config.get

        # Getting gas properties
        s.a_gas_dens = float(get('Properties_gas', 'a_gas_dens'))
        s.b_gas_dens = float(get('Properties_gas', 'b_gas_dens'))
        s.gas_visc = float(get('Properties_gas', 'gas_visc'))

        # Getting liquid properties
        s.liq_dens = float(get('Properties_liquid', 'liq_dens'))
        s.liq_visc = float(get('Properties_liquid', 'liq_visc'))

    def get_props_array(s):
        props_list = list()
        props_list.append(s.a_gas_dens)
        props_list.append(s.b_gas_dens)
        props_list.append(s.gas_visc)
        props_list.append(s.liq_dens)
        props_list.append(s.liq_visc)
        return np.array(props_list, dtype=float)

    def __str__(s):
        out_str = 'a_gas_dens ' + str(s.a_gas_dens)
        out_str += '\nb_gas_dens ' + str(s.b_gas_dens)
        out_str += '\ngas_visc ' + str(s.gas_visc)
        out_str += '\nliq_dens ' + str(s.liq_dens)
        out_str += '\nliq_visc ' + str(s.liq_visc)
        return out_str


if __name__ == '__main__':
    props = Props(config_file=sys.argv[1])
    print(props.get_props_array())
    print(props)




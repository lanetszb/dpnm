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

        # Getting gas properties PNM
        s.a_gas_dens = float(get('Properties_gas', 'a_gas_dens'))
        s.b_gas_dens = float(get('Properties_gas', 'b_gas_dens'))
        s.gas_visc = float(get('Properties_gas', 'gas_visc'))

        # Getting liquid properties PNM
        s.liq_dens = float(get('Properties_liquid', 'liq_dens'))
        s.liq_visc = float(get('Properties_liquid', 'liq_visc'))

        # Getting simulation properties PNM
        s.pressure_in = float(get('Properties_simulation', 'pressure_in'))
        s.pressure_out = float(get('Properties_simulation', 'pressure_out'))
        s.it_accuracy_PNM = float(get('Properties_simulation', 'it_accuracy'))

        # Getting diffusion properties
        s.time = float(get('Properties_diffusion', 'time'))
        s.time_step = float(get('Properties_diffusion', 'time_step'))
        s.length = float(get('Properties_diffusion', 'length'))
        s.radius = float(get('Properties_diffusion', 'radius'))
        s.eff_radius = float(get('Properties_diffusion', 'effective_radius'))
        s.grid_block_n = float(get('Properties_diffusion', 'grid_block_n'))
        s.conc_ini = float(get('Properties_diffusion', 'conc_ini'))
        s.diffusivity = float(get('Properties_diffusion', 'diffusivity'))
        s.it_accuracy = float(get('Properties_diffusion', 'iterative_accuracy'))

        # Getting langmuir coeffs
        s.langm_coeff = str(get('Langmuir_isotherm', 'langm_coeff'))
        s.langm_coeff = np.loadtxt(s.langm_coeff)

        # Getting matrix volume
        s.matrix_volume = float(get('Properties_matrix', 'matrix_volume'))

    def get_props_array(s):
        props_list = list()
        props_list.append(s.a_gas_dens)
        props_list.append(s.b_gas_dens)
        props_list.append(s.gas_visc)
        props_list.append(s.liq_dens)
        props_list.append(s.liq_visc)
        props_list.append(s.pressure_in)
        props_list.append(s.pressure_out)
        props_list.append(s.it_accuracy_PNM)
        return np.array(props_list, dtype=float)

    def get_diff_props_array(s):
        props_list = list()
        props_list.append(s.time)
        props_list.append(s.time_step)
        props_list.append(s.length)
        props_list.append(s.radius)
        props_list.append(s.eff_radius)
        props_list.append(s.grid_block_n)
        props_list.append(s.conc_ini)
        props_list.append(s.diffusivity)
        props_list.append(s.it_accuracy)
        return np.array(props_list, dtype=float)

    def __str__(s):
        out_str = 'a_gas_dens ' + str(s.a_gas_dens)
        out_str += '\nb_gas_dens ' + str(s.b_gas_dens)
        out_str += '\ngas_visc ' + str(s.gas_visc)
        out_str += '\nliq_dens ' + str(s.liq_dens)
        out_str += '\nliq_visc ' + str(s.liq_visc)
        out_str += '\npress_in ' + str(s.pressure_in)
        out_str += '\npress_out ' + str(s.pressure_out)
        out_str += '\nit_accuracy ' + str(s.it_accuracy_PNM)
        return out_str


if __name__ == '__main__':
    props = Props(config_file=sys.argv[1])
    print('PNM properties vector: ', props.get_props_array())
    print('\n')
    print('Diffusion properties vector: ', props.get_diff_props_array())
    print('\n')
    print(props)

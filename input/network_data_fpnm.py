import sys
import os
import pandas as pd
import configparser

import numpy as np
from input import plot_network_conns, plot_network_stats

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))


class Network_Data_Fpnm:
    def __init__(s, config_file):
        s.config = configparser.ConfigParser()
        s.config.read(config_file)

        # Read the path to text files with pores and throats

        s.fractures_data = str(s.config.get('PNData', 'fractures_data'))
        s.pores_data = str(s.config.get('PNData', 'pores_data'))
        s.boundary_pores = str(s.config.get('PNData', 'boundary_pores'))

        # Declare some future variables
        # Throat data
        s.fractures = None

        s.fractures_list = None
        s.fractures_widths = None
        s.fractures_lengths = None
        s.fractures_heights = None

        s.fracs_conn_ind_in = None
        s.fracs_conn_ind_out = None

        # Pore data
        s.pores = None
        s.pores_coords = None
        s.pores_list = None
        s.pores_radii = None

        s.inlet_pores = None
        s.outlet_pores = None

        s.pores_coords_x = None
        s.pores_coords_y = None
        s.pores_coords_z = None
        # boundary pores

        s.pores_left_x = None
        s.pores_right_x = None

        s.pores_front_y = None
        s.pores_back_y = None

        s.pores_top_z = None
        s.pores_bot_z = None

        # hydraulic conductance

        s.hydraulic_cond_coeff = None

    def process_fractures(s):
        s.fractures = pd.read_csv(s.fractures_data)
        s.fractures = s.fractures.sort_values(by=['pore_i', 'pore_j'])

        s.fracs_conn_ind_in = list(s.fractures['pore_i'])
        s.fracs_conn_ind_out = list(s.fractures['pore_j'])

        s.fractures_list = list(s.fractures['fracs_list'])
        s.fractures_heights = list(s.fractures['fracs_heights'])
        s.fractures_widths = list(s.fractures['fracs_widths'])
        s.fractures_lengths = list(s.fractures['fracs_lengths'])

        # hydraulic coeffs
        s.hydraulic_cond_coeff = list(s.fractures['hydr_cond'])

    def process_pores(s):
        s.pores = pd.read_csv(s.pores_data)
        s.boundary_pores = pd.read_csv(s.boundary_pores)

        s.pores_coords_x = list(s.pores['x_coord'])
        s.pores_coords_y = list(s.pores['y_coord'])
        s.pores_coords_z = list(s.pores['z_coord'])

        # s.pores_numbers = len(s.pores_coords_x)
        s.pores_list = list(s.pores['pores_list'])
        s.pores_radii = list(s.pores['pores_diameter'] / 2)

        # boundary pores
        s.pores_left_x = s.boundary_pores['pore_left_x']
        s.pores_right_x = s.boundary_pores['pore_right_x']

        s.pores_front_y = s.boundary_pores['pore_front_y']
        s.pores_back_y = s.boundary_pores['pore_back_y']

        s.pores_top_z = s.boundary_pores['pore_top_z']
        s.pores_bot_z = s.boundary_pores['pore_bot_z']


if __name__ == '__main__':
    network_data_cpp = Network_Data_Fpnm(config_file=sys.argv[1])

    network_data_cpp.process_fractures()
    network_data_cpp.process_pores()

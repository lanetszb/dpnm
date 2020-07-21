import sys
import os
import pandas as pd
import configparser

import numpy as np
from input import plot_network_conns, plot_network_stats

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))


class Network_Data_Conv:
    def __init__(s, config_file):
        s.config = configparser.ConfigParser()
        s.config.read(config_file)

        # Read the path to text files with pores and throats

        s.pore_throats = str(s.config.get('PNData', 'pore_throats'))
        s.pores_data = str(s.config.get('PNData', 'pores_data'))
        s.boundary_pores = str(s.config.get('PNData', 'boundary_pores'))

        # Declare some future variables
        # Throat data
        s.throats = None

        s.conn_ind = None
        s.throat_radius = None
        s.throat_length = None
        s.throat_number = None
        s.throat_list = None

        # Pore data
        s.pores = None
        s.pore_coords = None
        s.pore_number = None
        s.pore_list = None
        s.pore_radius = None

        s.inlet_pores = None
        s.outlet_pores = None
        # s.boundary_pores = None

        # boundary pores

        s.pore_left_x = None
        s.pore_right_x = None

        s.pore_front_y = None
        s.pore_back_y = None

        s.pore_top_z = None
        s.pore_bot_z = None

        # hydraulic conductance

        s.hydraulic_cond_coeff = None

    def process_throats(s):
        s.throats = pd.read_csv(s.pore_throats, index_col=0)
        s.throats = s.throats.sort_values(by=['pore_i', 'pore_j'])

        s.conn_ind_in = list(s.throats['pore_i'])
        s.conn_ind_out = list(s.throats['pore_j'])
        s.conn_ind = [[a, b] for a, b in zip(s.throats['pore_i'],
                                             s.throats['pore_j'])]

        s.throat_number = len(s.conn_ind_in)
        s.throat_list = np.arange(s.throat_number).tolist()

        s.throat_radius = list(s.throats['throat_height'] / 2)
        s.throat_width = list(s.throats['throat_width'])
        s.throat_length = list(s.throats['throat_length'])

        # hydraulic coeffs
        s.hydraulic_cond_coeff = list(s.throats['hydr_cond'])

    def process_pores(s):
        s.pores = pd.read_csv(s.pores_data, index_col=0)
        s.boundary_pores = pd.read_csv(s.boundary_pores)

        s.pore_coords_x = list(s.pores['x_coord'])
        s.pore_coords_y = list(s.pores['y_coord'])
        s.pore_coords_z = list(s.pores['z_coord'])

        s.pore_number = len(s.pore_coords_x)
        s.pore_list = np.arange(s.pore_number).tolist()
        s.pore_radius = list(s.pores['pore_diameter'] / 2)

        # boundary pores
        s.pore_left_x = s.boundary_pores['pore_left_x']
        s.pore_right_x = s.boundary_pores['pore_right_x']

        s.pore_front_y = s.boundary_pores['pore_front_y']
        s.pore_back_y = s.boundary_pores['pore_back_y']

        s.pore_top_z = s.boundary_pores['pore_top_z']
        s.pore_bot_z = s.boundary_pores['pore_bot_z']


if __name__ == '__main__':
    network_data_cpp = Network_Data_Conv(config_file=sys.argv[1])

    network_data_cpp.process_throats()
    network_data_cpp.process_pores()

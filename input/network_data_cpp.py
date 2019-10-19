import sys
import os
import pandas as pd
import configparser

import numpy as np
from input import plot_network_conns, plot_network_stats

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))


class Network_Data_Cpp:
    def __init__(s, config_file):
        s.config = configparser.ConfigParser()
        s.config.read(config_file)

        # Read the path to text files with pores and throats

        s.pore_throats = str(s.config.get('PNData', 'pore_throats'))
        s.pores_data = str(s.config.get('PNData', 'pores_data'))

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
        s.pore_conns = None
        s.pore_radius = None

        s.inlet_pores = None
        s.outlet_pores = None
        s.boundary_pores = None
        s.pore_per_row = None

        s.all_conns = None

    def process_throats(s):
        s.throats = pd.read_csv(s.pore_throats, index_col=0)
        s.throats = s.throats.sort_values(by=['pore_i', 'pore_j'])

        s.conn_ind_in = list(s.throats['pore_i'])
        s.conn_ind_out = list(s.throats['pore_j'])

        s.throat_number = len(s.conn_ind_in)
        s.throat_list = np.arange(s.throat_number).tolist()

        s.throat_radius = list(s.throats['throat_diameter'] / 2)
        s.throat_length = list(s.throats['throat_length'])

    def process_pores(s):
        s.pores = pd.read_csv(s.pores_data, index_col=0)

        s.pore_coords_x = list(s.pores['x_coord'])
        s.pore_coords_y = list(s.pores['y_coord'])
        s.pore_coords_z = list(s.pores['z_coord'])

        s.pore_number = len(s.pore_coords_x)
        s.pore_list = np.arange(s.pore_number).tolist()
        s.pore_radius = list(s.pores['pore_diameter'] / 2)

        s.conn_number = s.pores['conn_number']
        s.conn_number = list(s.conn_number)

    def process_pore_conns(s):

        pore_conns_row = s.pores['conn_indices']
        s.pore_conns = []

        for i in range(len(pore_conns_row)):
            pore_conns_row[i] = pore_conns_row[i].split(',')

        for i in range(len(s.pore_list)):
            for j in range(s.conn_number[i]):
                pore_conns_row[i][j] = int(pore_conns_row[i][j])

        for i in range(len(pore_conns_row)):
            pore_conns_row[i].sort()

        for i in range(len(pore_conns_row)):
            for j in range(s.conn_number[i]):
                s.pore_conns.append(int(pore_conns_row[i][j]))

    def process_pore_per_row(s):

        s.pore_per_row = s.pores['conn_indices']

        for i in range(len(s.pore_list)):
            for j in range(s.conn_number[i]):
                s.pore_per_row[i][j] = int(s.pore_per_row[i][j])

        s.pore_per_row = [[a] + b for a, b in zip(s.pore_list,
                                                     s.pore_per_row)]
        s.pore_per_row = np.array(
            [elem for single_list in s.pore_per_row for elem in
             single_list])

        s.pore_per_row = [int(i) for i in s.pore_per_row]


if __name__ == '__main__':
    network_data_cpp = Network_Data_Cpp(config_file=sys.argv[1])

    network_data_cpp.process_throats()
    network_data_cpp.process_pores()
    network_data_cpp.process_pore_conns()
    network_data_cpp.process_pore_per_row()
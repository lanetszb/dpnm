import sys
import os
import pandas as pd
import configparser

import numpy as np
from input import plot_network_conns, plot_network_stats

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))


class Network_Data:
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

        s.all_conns = None

    # Process throat data, create a list of connected pores

    def process_throats(s):
        s.throats = pd.read_csv(s.pore_throats, index_col=0)
        s.throats = s.throats.sort_values(by=['pore_i'])
        s.conn_ind = [[a, b] for a, b in zip(s.throats['pore_i'],
                                             s.throats['pore_j'])]

        s.throat_radius = s.throats['throat_height'] / 2
        s.throat_length = s.throats['throat_length']

    # Process pore data, create a list of pore coordinates

    def process_pores(s):
        s.pores = pd.read_csv(s.pores_data, index_col=0)
        s.pore_coords = [[a, b, c] for a, b, c in
                         zip(s.pores['x_coord'],
                             s.pores['y_coord'],
                             s.pores['z_coord'])]
        s.pore_number = len(s.pore_coords)
        s.pore_list = np.arange(len(s.pore_coords))
        s.pore_radius = s.pores['pore_diameter'] / 2

    # Identify boundary pores in any direction (input list of x, y or z coords)

    def get_boundary_pores(s, coord_list):

        s.boundary_pores = []

        x_min = min(coord_list)
        x_max = max(coord_list)

        s.inlet_pores = []
        s.outlet_pores = []

        for i in range(len(coord_list)):
            if coord_list[i] == x_min:
                s.inlet_pores.append(i)

        for i in range(len(coord_list)):
            if coord_list[i] == x_max:
                s.outlet_pores.append(i)

        s.boundary_pores = s.inlet_pores + s.outlet_pores

    # !Significantly slows down the whole pore processing
    def process_pore_conns(s):

        s.pore_conns = [[] for i in range(len(s.pore_list))]

        for i in range(len(s.pore_list)):
            for j in range(len(s.conn_ind)):
                if i == s.conn_ind[j][0]:
                    s.pore_conns[i].append(s.conn_ind[j][1])
                if i == s.conn_ind[j][1]:
                    s.pore_conns[i].append(s.conn_ind[j][0])
                if j == s.conn_ind:
                    break

    def get_all_conns(s):
        s.pore_list = s.pore_list.reshape(s.pore_number, 1)
        s.pore_list = s.pore_list.tolist()

        row_conns = []
        s.all_conns = []

        for i in range(s.pore_number):
            for j in range(len(s.pore_conns[i])):
                single_conn = list(s.pore_list[i] + [s.pore_conns[i][j]])
                single_conn.sort()
                row_conns.append(single_conn)
            s.all_conns.append(row_conns)
            row_conns = []

    # Print processed data in string format and print it to console

    def __str__(s):
        out_str = 'pore_throats ' + str(s.pore_throats) + '\n' + \
                  'pores_data ' + str(s.pores_data) + '\n'

        if s.throats is not None:
            out_str += '\nthroats ' + '\n' + str(s.throats) + '\n'
        if s.conn_ind is not None:
            out_str += '\nconn_ind ' + '\n' + str(s.conn_ind) + '\n'
        if s.throat_radius is not None:
            out_str += '\nthroat_radius ' + '\n' + str(s.throat_radius) + '\n'
        if s.throat_length is not None:
            out_str += '\nthroat_length ' + '\n' + str(s.throat_length) + '\n'

        if s.pores is not None:
            out_str += '\npores ' + '\n' + str(s.pores) + '\n'
        if s.pore_coords is not None:
            out_str += '\npore_coords ' + '\n' + str(s.pore_coords) + '\n'
        if s.pore_conns is not None:
            out_str += '\npore_conns ' + '\n' + str(s.pore_conns) + '\n'
        if s.pore_radius is not None:
            out_str += '\npore_radius ' + '\n' + str(s.pore_radius)
        if s.all_conns is not None:
            out_str += '\nall_conns ' + '\n' + str(s.all_conns)

        return out_str


if __name__ == '__main__':
    network_data = Network_Data(config_file=sys.argv[1])

    network_data.process_throats()
    network_data.process_pores()

    network_data.process_pore_conns()
    network_data.get_all_conns()

    print(network_data)

    conn_ind = network_data.conn_ind
    pore_coords = network_data.pore_coords
    pore_list = network_data.pore_list

    plot_network_conns(conn_ind, pore_coords, pore_list)

    property_to_plt = network_data.throat_radius
    plot_network_stats(property_to_plt)

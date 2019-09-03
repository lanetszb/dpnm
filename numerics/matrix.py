import sys
import os
import configparser
import numpy as np
import math
import pandas as pd

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from input import Props, Network_Data


class Matrix:
    def __init__(s, config_file):
        s.__config = configparser.ConfigParser()
        s.__config.read(config_file)

        s.props = Props(config_file)
        s.liq_visc = s.props.liq_visc

        s.network_stat = Network_Data(config_file)
        s.network_stat.process_throats()
        s.network_stat.process_pores()

        s.throat_radius = s.network_stat.throats['throat_diameter'] / 2
        s.conn_ind = np.array(s.network_stat.conn_ind)
        # s.conn_ind.sort()
        s.conn_number = s.network_stat.pores['conn_number']
        s.pore_coords = np.array(s.network_stat.pore_coords)

        s.pore_number = len(s.pore_coords)
        s.pore_list = np.arange(len(s.pore_coords))
        s.conn_indices = np.array(s.network_stat.conn_indices)

        s.throat_length = float
        s.matrix_coeff = type(s.throat_radius)
        s.matrix_row = []

    def center_distance(s):
        def xyz_distance(coord1, coord2):
            return np.sqrt(np.sum((coord1 - coord2) ** 2, axis=1))

        # Calculate center to center distances between the connected pores,
        # results have the same indexation as throats.
        #
        # Parameters:
        # pn = openpnm pore network object. Pore coordinates and
        # throat connections have to be known.

        p1 = s.pore_coords[s.conn_ind[:, 0]]
        p2 = s.pore_coords[s.conn_ind[:, 1]]
        return xyz_distance(p1, p2)

    def get_matrix_coeff(s):
        s.throat_length = matrix.center_distance()
        s.matrix_coeff = math.pi * matrix.throat_radius ** 4 / \
                         (8 * matrix.liq_visc * throat_length)

    #    def matrix_col(s):

    def get_matrix_row(s):

        for i in range(s.pore_number):
            for j in range(s.conn_number[i] + 1):
                s.matrix_row.append(i)

        s.matrix_row = np.array(s.matrix_row).flatten()

    def get_matrix_col(s):
        s.conn_indices = list(s.conn_indices)

        for i in range(len(s.conn_indices)):
            s.conn_indices[i] = s.conn_indices[i].split(',')
        for i in range(len(s.conn_indices)):
            for j in range(len(s.conn_indices[i])):
                s.conn_indices[i][j] = int(s.conn_indices[i][j])

        s.pore_list = s.pore_list.reshape(s.pore_number, 1)
        s.pore_list = s.pore_list.tolist()

        s.matrix_col = [a + b for a, b in zip(matrix.pore_list,
                                              matrix.conn_indices)]

        s.matrix_col = [y for x in s.matrix_col for y in x]

        s.matrix_col = np.array(s.matrix_col)

    def get_matrix_val(s):
        conn_ind = sorted(matrix.conn_ind, key=lambda a_entry: a_entry[1])
        conn_ind = np.array(conn_ind).tolist()

        temp = []
        temp1 = []

        for i in range(s.pore_number):
            for j in range(len(s.conn_indices[i])):
                temp.append(s.pore_list[i] + [s.conn_indices[i][j]])
            temp1.append(temp)
            temp = []

        print(temp1)

        temp1 = []



        # for i in range(9):
        #     for j in range(len(matrix.conn_indices[i])):
        #         x.append(matrix.pore_list[i] + [matrix.conn_indices[i][j]])


if __name__ == '__main__':
    matrix = Matrix(config_file=sys.argv[1])
    throat_length = matrix.center_distance()

    matrix.get_matrix_coeff()
    matrix.get_matrix_col()
    matrix.get_matrix_row()

    # matrix_coeff = math.pi * matrix.throat_radius ** 4 / \
    #                (8 * matrix.liq_visc * throat_length)

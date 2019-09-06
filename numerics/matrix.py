import sys
import os
import configparser
import numpy as np
import math
import pandas as pd

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from input import Props, Network_Data

from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve


class Matrix:
    def __init__(s, config_file):
        s.__config = configparser.ConfigParser()
        s.__config.read(config_file)

        s.props = Props(config_file)
        s.liq_visc = s.props.liq_visc

        s.network_stat = Network_Data(config_file)
        s.network_stat.process_throats()
        s.network_stat.process_pores()
        s.network_stat.boundary_pores()

        s.throat_radius = s.network_stat.throats['throat_diameter'] / 2
        s.conn_ind = np.array(s.network_stat.conn_ind)
        # s.conn_ind.sort()
        s.conn_number = s.network_stat.pores['conn_number']
        s.pore_coords = np.array(s.network_stat.pore_coords)

        s.pore_number = len(s.pore_coords)
        s.pore_list = np.arange(len(s.pore_coords))
        s.conn_indices = np.array(s.network_stat.pores['conn_indices'])

        s.throat_length = float
        s.matrix_coeff = type(s.throat_radius)
        s.matrix_row = []

        s.boundaries = s.network_stat.front_boundaries

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

    def get_matrix_row(s):

        for i in range(s.pore_number):
            x = []
            for j in range(s.conn_number[i] + 1):
                x.append(i)
            s.matrix_row.append(x)

        for i in range(len(s.boundaries)):
            for j in range(len(s.matrix_row)):
                if s.boundaries[i] == j:
                    s.matrix_row[j] = []
                    s.matrix_row[j].append(s.boundaries[i])

        s.matrix_row = [y for x in s.matrix_row for y in x]
        s.matrix_row = np.array(s.matrix_row)

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

        for i in range(len(s.boundaries)):
            for j in range(len(s.matrix_col)):
                if s.boundaries[i] == j:
                    s.matrix_col[j] = []
                    s.matrix_col[j].append(s.boundaries[i])

        print(s.matrix_col)

        s.matrix_col = [y for x in s.matrix_col for y in x]

        s.matrix_col = np.array(s.matrix_col)

    def get_matrix_val(s):
        conn_ind = np.array(s.conn_ind).tolist()
        print(conn_ind)

        temp = []
        s.temp1 = []

        for i in range(s.pore_number):
            for j in range(len(s.conn_indices[i])):
                x = list(s.pore_list[i] + [s.conn_indices[i][j]])
                x.sort()
                temp.append(x)
            s.temp1.append(temp)
            temp = []

        matrix_coeff = np.array(s.matrix_coeff)
        matrix_coeff = matrix_coeff.reshape(24, 1)
        matrix_coeff = matrix_coeff.tolist()

        s.conn_coeff = [[a, b] for a, b in
                        zip(conn_ind, matrix_coeff)]

        s.matrix_val = []
        temp = []

        for i in range(len(s.temp1)):
            for j in range(len(s.temp1[i])):
                for p in range(len(conn_ind)):
                    if s.temp1[i][j] == s.conn_coeff[p][0]:
                        temp.append(-1 * s.conn_coeff[p][1][0])
            s.matrix_val.append(list(temp))
            temp = []

        s.matrix_central = []
        for i in range(len(s.matrix_val)):
            s.matrix_central.append(-1 * sum(s.matrix_val[i]))

        for i in range(s.pore_number):
            s.matrix_val[i].insert(0, s.matrix_central[i])

        # s.matrix_val = [a + b for a, b in zip(s.matrix_central, s.matrix_val)]

        # x = [a + b for a, b in zip(s.pore_list,
        #                            s.conn_indices)]
        # for i in range(matrix.pore_number):
        #     x[i].sort()
        #
        # s.pore_index = []
        # for i in range(s.pore_number):
        #     for j in range(len(x[i])):
        #         if s.pore_list[i][0] == x[i][j]:
        #             s.pore_index.append(j)
        #
        # for i in range(s.pore_number):
        #     s.matrix_val[i].insert(s.pore_index[i], s.matrix_central[i])

        for i in range(len(s.boundaries)):
            for j in range(len(s.matrix_val)):
                if s.boundaries[i] == j:
                    s.matrix_val[j] = []
                    s.matrix_val[j].append(1)

        s.matrix_val = [y for x in s.matrix_val for y in x]
        s.matrix_val = np.array(matrix.matrix_val)


if __name__ == '__main__':
    matrix = Matrix(config_file=sys.argv[1])
    throat_length = matrix.center_distance()

    matrix.get_matrix_coeff()
    matrix.get_matrix_col()
    matrix.get_matrix_row()
    matrix.get_matrix_val()

    A = csr_matrix((matrix.matrix_val, (matrix.matrix_row, matrix.matrix_col)),
                   shape=(matrix.pore_number, matrix.pore_number)).toarray()
    B = np.array(
        [202650, 202650, 202650, 202650, 0, 0, 0, 0, 0, 0, 0, 0, 101325, 101325, 101325,
         101325])
    x = spsolve(A, B)

    # matrix_coeff = math.pi * matrix.throat_radius ** 4 / \
    #                (8 * matrix.liq_visc * throat_length)

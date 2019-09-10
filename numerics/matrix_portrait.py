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


class Matrix_Portrait:
    def __init__(s, config_file):
        s.__config = configparser.ConfigParser()
        s.__config.read(config_file)

        # Getting properties from corresponding class

        s.props = Props(config_file)
        s.liq_visc = s.props.liq_visc

        # Getting network data from corresponding class

        s.network_data = Network_Data(config_file)
        s.network_data.process_throats()
        s.network_data.process_pores()

        s.throat_radius = s.network_data.throat_radius
        s.throat_length = s.network_data.throat_length

        s.conn_ind = s.network_data.conn_ind

        s.conn_number = s.network_data.pores['conn_number']
        s.pore_coords = s.network_data.pore_coords

        s.pore_number = s.network_data.pore_number
        s.pore_list = s.network_data.pore_list

        s.front_pores = s.network_data.pores['x_coord']

        s.network_data.get_boundary_pores(s.front_pores)
        s.boundary_pores = s.network_data.boundary_pores

        s.matrix_coeff = None
        s.matrix_row = None
        s.matrix_col = None
        s.vector_B = None

        s.network_data.process_pore_conns()
        s.network_data.get_all_conns()

        s.pore_conns = s.network_data.pore_conns
        s.all_conns = s.network_data.all_conns

    def get_matrix_coeff(s):

        s.matrix_coeff = math.pi * s.throat_radius ** 4 / \
                         (8 * s.liq_visc * s.throat_length)

    def get_matrix_row(s):

        s.matrix_row = []

        for i in range(s.pore_number):
            pores_per_row = []
            for j in range(s.conn_number[i] + 1):
                pores_per_row.append(i)
            s.matrix_row.append(pores_per_row)

        for idx, value in enumerate(s.boundary_pores):
            s.matrix_row[value] = [value]

        s.matrix_row = [y for x in s.matrix_row for y in x]
        s.matrix_row = np.array(s.matrix_row)

    def get_matrix_col(s):

        s.pore_list = s.pore_list.reshape(s.pore_number, 1)
        s.pore_list = s.pore_list.tolist()

        s.matrix_col = [a + b for a, b in zip(s.pore_list,
                                              s.pore_conns)]

        for idx, value in enumerate(s.boundary_pores):
            s.matrix_col[value] = [value]

        s.matrix_col = [y for x in s.matrix_col for y in x]
        s.matrix_col = np.array(s.matrix_col)

    def get_matrix_val(s):
        conn_ind = np.array(s.conn_ind).tolist()

        matrix_coeff = np.array(s.matrix_coeff)
        matrix_coeff = matrix_coeff.reshape(len(s.conn_ind), 1)
        matrix_coeff = matrix_coeff.tolist()

        conn_coeff = [[a, b] for a, b in zip(conn_ind, matrix_coeff)]

        s.matrix_val = []
        coeff_in_row = []

        for i in range(len(s.all_conns)):
            for j in range(len(s.all_conns[i])):
                for k in range(len(conn_ind)):
                    if s.all_conns[i][j] == conn_coeff[k][0]:
                        coeff_in_row.append(-1 * conn_coeff[k][1][0])
            s.matrix_val.append(list(coeff_in_row))
            coeff_in_row = []

        s.matrix_central = []
        for i in range(len(s.matrix_val)):
            s.matrix_central.append(-1 * sum(s.matrix_val[i]))

        for i in range(s.pore_number):
            s.matrix_val[i].insert(0, s.matrix_central[i])

        for idx, value in enumerate(s.boundary_pores):
            s.matrix_val[value] = [1]

        s.matrix_val = [y for x in s.matrix_val for y in x]
        s.matrix_val = np.array(s.matrix_val)

    def get_vector_B(s, P_in, P_out):

        s.vector_B = list([0] * s.pore_number)

        for idx, value in enumerate(s.network_data.inlet_pores):
            s.vector_B[value] = P_in

        for idx, value in enumerate(s.network_data.outlet_pores):
            s.vector_B[value] = P_out

if __name__ == '__main__':
    matrix_portrait = Matrix_Portrait(config_file=sys.argv[1])

    matrix_portrait.get_matrix_coeff()
    matrix_portrait.get_matrix_row()
    matrix_portrait.get_matrix_col()
    matrix_portrait.get_matrix_val()

    P_in = 202650
    P_out = 101325

    matrix_portrait.get_vector_B(202650, 101325)
    #
    A = csr_matrix((matrix_portrait.matrix_val, (matrix_portrait.matrix_row, matrix_portrait.matrix_col)),
                   shape=(matrix_portrait.pore_number, matrix_portrait.pore_number)).toarray()

    x = spsolve(A, matrix_portrait.vector_B)

    # matrix_coeff = math.pi * matrix.throat_radius ** 4 / \
    #                (8 * matrix.liq_visc * throat_length)

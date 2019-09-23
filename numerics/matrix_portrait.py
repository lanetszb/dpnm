import sys
import os
import configparser
import numpy as np
import math
import copy

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

        # Getting network data from corresponding class

        s.network_data = Network_Data(config_file)

        s.network_data.process_throats()
        s.network_data.process_pores()
        s.network_data.process_pore_conns()
        s.network_data.get_all_conns()

        s.network_data.get_boundary_pores(s.network_data.pores['x_coord'])

        s.matrix_coeff = None
        s.matrix_row = None
        s.matrix_col = None
        s.matrix_val = None
        s.vector_B = None

    def get_matrix_coeff(s):

        s.matrix_coeff = math.pi * s.network_data.throat_radius ** 4 / \
                         (8 * s.props.liq_visc * s.network_data.throat_length)

    def get_matrix_row(s):

        s.matrix_row = []

        for i in range(s.network_data.pore_number):
            pores_per_row = []
            for j in range(s.network_data.pores['conn_number'][i] + 1):
                pores_per_row.append(i)
            s.matrix_row.append(pores_per_row)

        for idx, value in enumerate(s.network_data.boundary_pores):
            s.matrix_row[value] = [value]

        s.matrix_row = [y for x in s.matrix_row for y in x]
        s.matrix_row = np.array(s.matrix_row)

    def get_matrix_col(s):

        s.network_data.pore_list = np.array(s.network_data.pore_list).tolist()

        s.matrix_col = [a + b for a, b in zip(s.network_data.pore_list,
                                              s.network_data.pore_conns)]

        for idx, value in enumerate(s.network_data.boundary_pores):
            s.matrix_col[value] = [value]

        s.matrix_col = [y for x in s.matrix_col for y in x]
        s.matrix_col = np.array(s.matrix_col)

    def get_matrix_val(s):
        conn_ind = np.array(s.network_data.conn_ind).tolist()

        matrix_coeff = np.array(s.matrix_coeff)
        matrix_coeff = matrix_coeff.reshape(len(conn_ind), 1)
        matrix_coeff = matrix_coeff.tolist()

        conn_coeff = [[a, b] for a, b in zip(conn_ind, matrix_coeff)]

        s.matrix_val = []
        coeff_in_row = []

        for i in range(len(s.network_data.all_conns)):
            for j in range(len(s.network_data.all_conns[i])):
                for k in range(len(conn_ind)):
                    if s.network_data.all_conns[i][j] == conn_coeff[k][0]:
                        coeff_in_row.append(-1 * conn_coeff[k][1][0])
            s.matrix_val.append(list(coeff_in_row))
            coeff_in_row = []


        matrix_central = []
        for i in range(len(s.matrix_val)):
            matrix_central.append(-1 * sum(s.matrix_val[i]))

        for i in range(s.network_data.pore_number):
            s.matrix_val[i].insert(0, matrix_central[i])

        for idx, value in enumerate(s.network_data.boundary_pores):
            s.matrix_val[value] = [1]

        s.matrix_val = [y for x in s.matrix_val for y in x]
        s.matrix_val = np.array(s.matrix_val)

    def get_vector_B(s, P_in, P_out):

        s.vector_B = list([0] * s.network_data.pore_number)

        for idx, value in enumerate(s.network_data.inlet_pores):
            s.vector_B[value] = P_in

        for idx, value in enumerate(s.network_data.outlet_pores):
            s.vector_B[value] = P_out

    def get_flow_rate(s, pressure_list):

        inlet_conns = []
        for i in range(len(s.network_data.conn_ind)):
            for j in range(len(s.network_data.inlet_pores)):
                if s.network_data.inlet_pores[j] == s.network_data.conn_ind[i][0]:
                    inlet_conns.append(matrix_portrait.network_data.conn_ind[i])

        dp = []

        for idx, value in enumerate(inlet_conns):
            dp.append(pressure[value[0]] - pressure[value[1]])

        flow_rates = s.matrix_coeff[0:len(inlet_conns)] * dp
        flow_rates = np.array(flow_rates)

        return np.sum(flow_rates)


if __name__ == '__main__':
    matrix_portrait = Matrix_Portrait(config_file=sys.argv[1])

    matrix_portrait.get_matrix_coeff()
    matrix_portrait.get_matrix_row()
    matrix_portrait.get_matrix_col()
    matrix_portrait.get_matrix_val()

    P_in = 202650
    P_out = 101325

    matrix_portrait.get_vector_B(202650, 101325)

    A = csr_matrix((matrix_portrait.matrix_val,
                    (matrix_portrait.matrix_row,
                     matrix_portrait.matrix_col)),
                   shape=(matrix_portrait.network_data.pore_number,
                          matrix_portrait.network_data.pore_number)).toarray()

    x = spsolve(A, matrix_portrait.vector_B)

    pressure = list(x)

    total_flow = matrix_portrait.get_flow_rate(pressure)

    x_min = min(matrix_portrait.network_data.pores['x_coord'])
    x_max = max(matrix_portrait.network_data.pores['x_coord'])
    x_length = x_max - x_min
    area = x_length ** 2

    k = total_flow * matrix_portrait.props.liq_visc * x_length / (
            area * (P_in - P_out))

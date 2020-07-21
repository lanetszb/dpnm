import sys
import os
import numpy as np
import math

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from input import Props, Network_Data_Python

from scipy.sparse import csr_matrix


class Matrix_Portrait:
    def __init__(s, props, network_data):
        s.props = props
        s.network_data = network_data

        s.matrix_coeff = None
        s.row = None
        s.col = None
        s.val = None
        s.A = None
        s.vector_B = None

    def get_matrix_coeff(s):

        liq_visc = s.props.liq_visc
        radius = s.network_data.throat_radius
        length = s.network_data.throat_length

        s.matrix_coeff = math.pi * radius ** 4 / 8 / liq_visc / length

    def get_matrix_row(s):

        s.row = []

        pore_n = s.network_data.pore_number
        conn_n = s.network_data.pores['conn_number']
        boundary_pores = s.network_data.boundary_pores

        for i in range(pore_n):
            pores_per_row = []
            for j in range(conn_n[i] + 1):
                pores_per_row.append(i)
            s.row.append(pores_per_row)

        for idx, value in enumerate(boundary_pores):
            s.row[value] = [value]

        s.row = [y for x in s.row for y in x]
        s.row = np.array(s.row)

    def get_matrix_col(s):

        pore_list = np.array(s.network_data.pore_list).tolist()
        pore_conns = s.network_data.pore_conns
        boundary_pores = s.network_data.boundary_pores

        s.col = [a + b for a, b in zip(pore_list, pore_conns)]

        for idx, value in enumerate(boundary_pores):
            s.col[value] = [value]

        s.col = [y for x in s.col for y in x]
        s.col = np.array(s.col)

    def get_matrix_val(s):
        conn_ind = np.array(s.network_data.conn_ind).tolist()
        all_conns = s.network_data.all_conns
        pore_n = s.network_data.pore_number


        matrix_coeff = np.array(s.matrix_coeff)
        matrix_coeff = matrix_coeff.reshape(len(conn_ind), 1)
        matrix_coeff = matrix_coeff.tolist()

        conn_coeff = [[a, b] for a, b in zip(conn_ind, matrix_coeff)]

        s.val = []
        coeff_in_row = []

        for i in range(len(all_conns)):
            for j in range(len(all_conns[i])):
                for k in range(len(conn_ind)):
                    if all_conns[i][j] == conn_coeff[k][0]:
                        coeff_in_row.append(-1 * conn_coeff[k][1][0])
            s.val.append(list(coeff_in_row))
            coeff_in_row = []

        matrix_central = []
        for i in range(len(s.val)):
            matrix_central.append(-1 * sum(s.val[i]))

        for i in range(pore_n):
            s.val[i].insert(0, matrix_central[i])

        for idx, value in enumerate(s.network_data.boundary_pores):
            s.val[value] = [1]

        s.val = [y for x in s.val for y in x]
        s.val = np.array(s.val)

    def get_matrix_portrait(s):

        pore_n = s.network_data.pore_number

        s.A = csr_matrix((s.val, (s.row, s.col)),
                         shape=(pore_n, pore_n)).toarray()


if __name__ == '__main__':
    props = Props(config_file=sys.argv[1])

    network_data = Network_Data_Python(config_file=sys.argv[1])
    network_data.process_throats()
    network_data.process_pores()
    network_data.process_pore_conns()
    network_data.get_all_conns()
    network_data.get_boundary_pores(network_data.pores['x_coord'])

    matrix_portrait = Matrix_Portrait(props, network_data)

    matrix_portrait.get_matrix_coeff()
    matrix_portrait.get_matrix_row()
    matrix_portrait.get_matrix_col()
    matrix_portrait.get_matrix_val()

    matrix_portrait.get_matrix_portrait()
    print(matrix_portrait.matrix_coeff)
    print(matrix_portrait.A)

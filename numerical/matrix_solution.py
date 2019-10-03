import sys
import os
import numpy as np


current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from scipy.sparse.linalg import spsolve
from time import perf_counter

from numerics import Matrix_Portrait
from input import Props, Network_Data


class Matrix_Solution:
    def __init__(s, matrix_portrait, network_data):
        s.matrix_portrait = matrix_portrait
        s.network_data = network_data

    def get_vector_B(s, pressure_in, pressure_out):
        pore_n = s.network_data.pore_number
        inlet_pores = s.network_data.inlet_pores
        outlet_pores = s.network_data.outlet_pores

        s.vector_B = list([0] * pore_n)

        for idx, value in enumerate(inlet_pores):
            s.vector_B[value] = pressure_in

        for idx, value in enumerate(outlet_pores):
            s.vector_B[value] = pressure_out

    def get_flow_rate(s, pressure_list):
        conn_ind = s.network_data.conn_ind
        inlet_pores = s.network_data.inlet_pores
        matrix_coeff = s.matrix_portrait.matrix_coeff

        inlet_conns = []
        for i in range(len(conn_ind)):
            for j in range(len(inlet_pores)):
                if inlet_pores[j] == conn_ind[i][0]:
                    inlet_conns.append(network_data.conn_ind[i])

        dp = []

        for idx, value in enumerate(inlet_conns):
            dp.append(pressure_list[value[0]] - pressure_list[value[1]])

        flow_rates = matrix_coeff[0:len(inlet_conns)] * dp
        flow_rates = np.array(flow_rates)

        return np.sum(flow_rates)


if __name__ == '__main__':

    t_start = perf_counter()

    props = Props(config_file=sys.argv[1])

    network_data = Network_Data(config_file=sys.argv[1])

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

    matrix_solution = Matrix_Solution(matrix_portrait, network_data)

    pressure_in = 300000
    pressure_out = 200000

    matrix_solution.get_vector_B(pressure_in, pressure_out)
    matrix_portrait.get_matrix_portrait()

    pressure = spsolve(matrix_portrait.A, matrix_solution.vector_B)
    total_flow = matrix_solution.get_flow_rate(pressure)
    #
    x_min = min(matrix_portrait.network_data.pores['x_coord'])
    x_max = max(matrix_portrait.network_data.pores['x_coord'])
    x_length = x_max - x_min
    area = x_length ** 2

    k = total_flow * matrix_portrait.props.liq_visc * x_length / (
            area * (pressure_in - pressure_out))

    t_stop = perf_counter()
    print("Elapsed time in seconds:", t_stop-t_start)

import sys
import os
import pandas as pd
import configparser

import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

sys.path.append(
    os.path.join(current_path, '../../../../../../projects/pmeal/OpenPNM/'))
sys.path.append(
    os.path.join(current_path, '../../../../../../projects/pmeal/porespy/'))

import scipy as sp
import openpnm as op


class Network_Data:
    def __init__(s, config_file):
        s.config = configparser.ConfigParser()
        s.config.read(config_file)

        s.pore_throats = str(s.config.get('PNData', 'pore_throats'))
        s.pores_data = str(s.config.get('PNData', 'pores_data'))

        def __str__(s):
            out_str = 'pore_throats ' + str(s.pore_throats) + '\n' + \
                      'pores_data ' + str(s.pores_data)
            return out_str

        s.throats = None
        s.conn_ind = None
        s.pores = None
        s.pore_coords = None
        s.conn_indices = None


    def process_throats(s):
        s.throats = pd.read_csv(s.pore_throats, index_col=0)
        s.conn_ind = [[a, b] for a, b in zip(s.throats['pore_i'],
                                             s.throats['pore_j'])]

    def process_pores(s):
        s.pores = pd.read_csv(s.pores_data, index_col=0)
        s.pore_coords = [[a, b, c] for a, b, c in
                         zip(s.pores['x_coord'],
                             s.pores['y_coord'],
                             s.pores['z_coord'])]


    def plot(s, conn_ind, pore_coords):
        pn1 = op.network.GenericNetwork(conn_ind, pore_coords)
        plt_connects = op.topotools.plot_connections(pn1)
        plt_coords = op.topotools.plot_coordinates(pn1, fig=plt_connects,
                                                   color='r', s=100)
        plt.show()

    def boundary_pores(s):

        s.front_boundaries = []

        x_min = min(s.pores['x_coord'])
        x_max = max(s.pores['x_coord'])

        s.inlet_pores = []
        s.outlet_pores = []

        for i in range(len(s.pores['x_coord'])):
            if s.pores['x_coord'][i] == x_min:
                s.inlet_pores.append(i)

        for i in range(len(s.pores['x_coord'])):
            if s.pores['x_coord'][i] == x_max:
                s.outlet_pores.append(i)

        s.front_boundaries = s.inlet_pores + s.outlet_pores

    def __str__(s):
        out_str = super().__str__()
        if s.throats is not None:
            out_str += '\nthroats ' + '\n' + str(s.throats)
        if s.conn_ind is not None:
            out_str += '\nconn_ind ' + '\n' + str(s.conn_ind)
        if s.pores is not None:
            out_str += '\npores ' + '\n' + str(s.pores)
        if s.pore_coords is not None:
            out_str += '\npore_coords ' + '\n' + str(s.pore_coords)
        return out_str


if __name__ == '__main__':
    network_data = Network_Data(config_file=sys.argv[1])
    network_data.process_throats()
    network_data.process_pores()

    conns = network_data.conn_ind
    coords = network_data.pore_coords
    network_data.plot(conns, coords)

    print(network_data)

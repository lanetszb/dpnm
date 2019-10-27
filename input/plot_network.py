import sys
import os

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

sys.path.append(
    os.path.join(current_path, '../../../../../../tmp/pmeal/OpenPNM/'))
sys.path.append(
    os.path.join(current_path, '../../../../../../tmp/pmeal/porespy/'))

import openpnm as op


def plot_network_conns(conn_ind, pore_coords, pore_list):

    fig = plt.figure()
    ax = Axes3D(fig)

    pore_network = op.network.GenericNetwork(conn_ind, pore_coords)
    plt_connects = op.topotools.plot_connections(pore_network, fig=fig)
    plt_coords = op.topotools.plot_coordinates(pore_network,
                                               fig=plt_connects,
                                               color='r', s=150)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    pore_list = [y for x in pore_list for y in x]

    for coord in range(len(pore_coords)):
        ax.text(pore_coords[coord][0], pore_coords[coord][1],
                pore_coords[coord][2], pore_list[coord])

    fig.show()


def plot_network_stats(property_list):

    fig2 = plt.figure()
    # Plot Histogram on x
    x = property_list
    plt.hist(x, bins=50, ec='black')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.gca().set(title='Throat radius', ylabel='Frequency')

    fig2.autofmt_xdate()

    fig2.show()


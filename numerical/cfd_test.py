import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from numerical import PropsCpp
from numerical import LocalCpp
from numerical import ConvectiveCpp
from numerical import EquationCpp

from numerical import PropsPNMCpp
from numerical import NetworkDataCpp
from numerical import EquationPNM
from numerical import DiffusionPNM

from input import Props
from input import Network_Data_Cpp

props = Props(config_file=sys.argv[1])
#
props_array = props.get_diff_props_array()

props_diff = props.get_props_array()

props_cpp = PropsCpp(props_array, props.langm_coeff)
print(props_cpp)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++

# local_cpp = LocalCpp(props_array, props.langm_coeff)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++

# convective_cpp = ConvectiveCpp(props_array, props.langm_coeff)
#
# equation_cpp = EquationCpp(props_array, props.langm_coeff)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# radius = 4.86487966261438e-06
# eff_radius = radius * 5
# thr_length = 0.9999723657e-4
# conc_wall = 50
#
# local_cpp.calculate_alpha(props_cpp.time, radius, eff_radius, thr_length)
#
# equation_cpp.cfdProcedure(conc_wall, radius, eff_radius, thr_length)
#
# conc = equation_cpp.getConc()
# radius_curr = local_cpp.radius_curr
#
# flow_rate = equation_cpp.getFlowRate()
#
# grid_centers = []
# for i in range(len(radius_curr)-1):
#     grid_centers.append((radius_curr[i+1] +
#                          radius_curr[i])/2)
#
# plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
# plt.gca().set(title='Concentration distribution, t = 0.1 (sec)',
#               xlabel='Radius (m)',
#               ylabel='Concentration (kg/m3)')
#
# print(conc)
# plt.plot(grid_centers, conc)
#
# plt.ylim(0, 50)
# plt.show()

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# print("\n")
# #
props_pnm = props.get_props_array()
props_pnm_cpp = PropsPNMCpp(props_pnm)
print(props_pnm_cpp)

network_data_cpp = Network_Data_Cpp(config_file=sys.argv[1])
network_data_cpp.process_throats()
network_data_cpp.process_pores()
network_data_cpp.process_pore_conns()
network_data_cpp.process_pore_per_row()

thrList = network_data_cpp.throat_list
tr = network_data_cpp.throat_radius
tl = network_data_cpp.throat_length

conn_in = network_data_cpp.conn_ind_in
conn_out = network_data_cpp.conn_ind_out

pc_x = network_data_cpp.pore_coords_x
pc_y = network_data_cpp.pore_coords_y
pc_z = network_data_cpp.pore_coords_z

pr = network_data_cpp.pore_radius

pl = network_data_cpp.pore_list
p_conn = network_data_cpp.pore_conns
conn_numb = network_data_cpp.conn_number
ppr = network_data_cpp.pore_per_row

plx = network_data_cpp.pore_left_x
prx = network_data_cpp.pore_right_x

nd_cpp = NetworkDataCpp(thrList, tr, tl, conn_in, conn_out, pc_x, pc_y, pc_z,
                        pr, pl, p_conn, conn_numb, ppr, plx, prx)

eq_pnm = EquationPNM(props_pnm, thrList, tr, tl, conn_in, conn_out, pc_x, pc_y,
                     pc_z, pr, pl, p_conn, conn_numb, ppr, plx, prx)
#
# diff_pnm = DiffusionPNM(props_pnm, thrList, tr, tl, conn_in,
#                         conn_out, pc_x, pc_y, pc_z, pr, pl, p_conn, conn_numb,
#                         ppr, props_array, props.langm_coeff)

# =============================================================================
# Figure 1 (Avg Pore Pressure and Avg Concentration)
# t = np.arange(0, props.time, props.time_step)
# data1 = diff_pnm.get_pressure_av()
# data2 = diff_pnm.get_conc_av()
#
# fig, ax1 = plt.subplots()
#
# color = 'tab:red'
# ax1.set_xlabel('time (s)', fontsize='large')
# ax1.set_ylabel('Average Pore Pressure (Pa)', color=color)
# ax1.plot(t, data1, color=color, label='Average Pore Pressure')
# ax1.tick_params(axis='y', labelcolor=color)
#
# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#
# color = 'tab:blue'
# ax2.set_ylabel('Avg Concentration (kg/m3)',
#                color=color, fontsize='large')  # we already handled the x-label with ax1
# ax2.plot(t, data2, color=color, label='Avg Matrix Concentration')
# ax2.tick_params(axis='y', labelcolor=color)
#
# fig.legend(loc="upper center", fontsize='x-large')
# fig.tight_layout()  # otherwise the right y-label is slightly clipped
# plt.show()

# =============================================================================
# Figure 2 (Total Flow Rate)
# t = np.arange(0, props.time, props.time_step)
# data1 = diff_pnm.get_pressure_av()
# data2 = diff_pnm.get_flow_pores_out()
# data3 = diff_pnm.get_flow_pores_in()
#
# fig, ax1 = plt.subplots()
#
# color = 'tab:red'
# ax1.set_xlabel('time (s)')
# ax1.set_ylabel('Average Pore Pressure (Pa)', color=color)
# ax1.set_ylim([250000, 300000])
#
# ax1.plot(t, data1, color=color, label='Average Pore Pressure')
# ax1.tick_params(axis='y', labelcolor=color)
#
# # ax1.legend(loc = 0)
#
# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#
# color = 'tab:blue'
# ax2.set_ylabel('Total Flow Rate (m/sec)',
#                color=color)  # we already handled the x-label with ax1
# ax2.plot(t, data2, color=color, label='Outlet Flow Rate')
# ax2.plot(t, data3, color='tab:green', label='Inlet Flow Rate')
# # ax2.plot(t, data2-data3, color=color)
#
# ax2.tick_params(axis='y', labelcolor=color)
#
# # ax2.legend(loc=0)
#
#
# fig.legend(loc="upper center", fontsize='x-large')
# fig.tight_layout()  # otherwise the right y-label is slightly clipped
# plt.show()

# Figure 3 (Langmuir isotherm and density)
# p_av = diff_pnm.get_pressure_av()
# p_av = np.arange(200000, 600000, 20000)
# # a_dens (kg/m3) is a coefficient for equation density = a * P + b
# a_dens = 6.71079e-06
# # b_dens (kg/m3/Pa) is b coefficient for equation density = a * P + b
# b_dens = -2.37253E-02
# rho_const = 1.653
#
#
# def calc_langm_conc(pressure):
#     langm_conc = 0
#     for i in range(len(props.langm_coeff)):
#         langm_conc += props.langm_coeff[i] * pressure ** i
#     return langm_conc
#
#
# langm_conc_list = []
# for i in range(len(p_av)):
#     langm_conc_list.append(calc_langm_conc(p_av[i]))
#
# density_list = []
# for i in range(len(p_av)):
#     density_list.append(a_dens * p_av[i] + b_dens)
#
# density_const_list = []
# for i in range(len(p_av)):
#     density_const_list.append(rho_const)
#
# fig, ax1 = plt.subplots()
#
# color = 'tab:red'
# ax1.set_xlabel('pressure (Pa)')
# ax1.set_ylabel('Langmuir Concentration (Pa)', color=color)
# ax1.set_ylim([0.7, 1.5])
#
# ax1.plot(p_av, langm_conc_list, color=color, label='Langmuir Concentration')
# ax1.tick_params(axis='y', labelcolor=color)
#
# # ax1.legend(loc = 0)
#
# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#
# color = 'tab:blue'
# ax2.set_ylabel('Gas Density (kg/m3)',
#                color=color)  # we already handled the x-label with ax1
# ax2.plot(p_av, density_const_list, color=color, label='Constant Density')
# ax2.plot(p_av, density_list, color='tab:green', label='Real Density')
# # ax2.plot(t, data2-data3, color=color)
#
# ax2.tick_params(axis='y', labelcolor=color)
#
# # ax2.legend(loc=0)
#
#
# fig.legend(loc="upper center", fontsize='x-large')
# fig.tight_layout()  # otherwise the right y-label is slightly clipped
# plt.show()
#
# # Figure 4 (Pressure and Concentration by Length)
# length_coords = list(dict.fromkeys(pc_x))
#
# pore_by_coord = [[] for x in range(len(length_coords))]
#
# for i in range(len(length_coords)):
#     for j in range(len(pc_x)):
#         if length_coords[i] == pc_x[j]:
#             pore_by_coord[i].append(j)
#
# for i in range(len(length_coords)):
#     length_coords[i] = [length_coords[i]]
#
# pores_by_coords = [[a] + [b] for a, b in zip(length_coords,
#                                            pore_by_coord)]
#
# pressure_by_pore = diff_pnm.get_pressure_by_pore()

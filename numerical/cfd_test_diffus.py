import sys
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from numerical import PropsDiffusion
from numerical import LocalDiffusion
from numerical import ConvectiveDiffusion
from numerical import EquationDiffusion

from numerical import PropsPNMCpp
from numerical import NetworkDataCpp
from numerical import EquationPNM
from numerical import DiffusionPNM

from input import Props
from input import Network_Data_Cpp

from output import plot_x_y
from output import plot_x_ymult
from matplotlib import rc

rc('text', usetex=True)

props = Props(config_file=sys.argv[1])

props_diff_vector = props.get_diff_props_array()
#
# props_diff = PropsDiffusion(props_diff_vector)
# print(props_diff, '\n')
#
# props_diff.print_props_vector()
# print('\n')

#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++

# convective_diff = ConvectiveDiffusion(props_diff_vector)
#
# local_diff = LocalDiffusion(props_diff_vector)
#
# equation_diff = EquationDiffusion(props_diff_vector)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# time_step = props_diff.time_step
# radius = props_diff.radius
# eff_radius = props_diff.eff_radius
# thr_length = props_diff.length
# grid_block_n = props_diff.grid_block_n
# diffusivity = props_diff.diffusivity
# conc_wall = 2.0
#
# #
# local_diff.calc_vol_cylinder(radius, eff_radius,
#                              grid_block_n, thr_length)
# #
# local_diff.calc_vol_cartesian(radius, eff_radius,
#                               thr_length, radius)
# #
# # # volume_list = local_diff.vol_cylindr
# volume_list = local_diff.vol_cartes
# #
# local_diff.calc_vol_cartesian(radius, eff_radius, thr_length, eff_radius / 10)
# #
# #
# # convective_diff.calc_omega_cylindr(thr_length)
# #
# convective_diff.calc_omega_cartes(radius, thr_length)
# #
# # # omega_list = convective_diff.omega_cylindr
# omega_list = convective_diff.omega_cartes
# #
# # #
# local_diff.calculate_alpha(time_step, volume_list)
# #
# convective_diff.calculate_beta(radius, eff_radius, thr_length,
#                                diffusivity, grid_block_n, omega_list)

# equation_diff.cfd_procedure_one_step(conc_wall, radius,
#                             eff_radius, thr_length,
#                             volume_list, omega_list)
#
# x = (eff_radius - radius) / (10 - 1)
#
# equation_diff.cfd_procedure(1, conc_wall, radius - x / 2,
#                             eff_radius + x / 2, thr_length,
#                             volume_list, omega_list)


# equation_diff.cfd_procedure(1, conc_wall, radius,
#                             eff_radius, thr_length,
#                             volume_list, omega_list)
# # #
# conc = equation_diff.getConc()
#
# radius_curr = local_diff.radius_curr
#
# flow_rate = equation_diff.getFlowRate()
# #
# grid_centers = []
# for i in range(len(radius_curr) - 1):
#     grid_centers.append((radius_curr[i + 1] +
#                          radius_curr[i]) / 2)
#
# # Diffusion Analytical Solution
# conc_analyt = []
#
# eff_radius_analyt = 2.43243983130719e-05
# radius_analyt = 4.86487966261438e-06
#
# # L = eff_radius - radius
# L = eff_radius_analyt - radius_analyt
#
# grid_block_n_analyt = 300
# dX = L / grid_block_n_analyt
#
# # dX = L / grid_block_n
# # dX = L / 100
#
# grid_centers_analyt = []
# for i in range(grid_block_n_analyt + 1):
#     grid_centers_analyt.append(radius_analyt + i * dX)
#     # grid_centers_analyt.append(radius + i * dX + dX / 2)
#
# conc_ini = 2.0
# conc_out = 10.0
# D = props_diff.diffusivity
# t = props_diff.time
# dt = props_diff.time_step
#
# for i in range(grid_block_n_analyt + 1):
#     conc_it = conc_out + (conc_ini - conc_out) * math.erf(
#         i * dX / 2 / math.sqrt(D * t))
#     conc_analyt.append(conc_it)
#
# conc_analyt.reverse()
#
# plot_x_y(grid_centers_analyt, conc_analyt, x_name='Radius (m)',
#          y_name='Concentration (kg/m3)',
#          graph_name='Concentration distribution',
#          line_type='-')

# ===================================================================

# plot_x_y(grid_centers, conc, x_name='Radius (m)',
#          y_name='Concentration (kg/m3)',
#          graph_name='Concentration distribution',
#          line_type='--',
#          lw=0.3, marker='s', markersize=4)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#
# print("\n")
#
props_pnm = props.get_props_array()
props_pnm_cpp = PropsPNMCpp(props_pnm)
# print(props_pnm_cpp)
#
network_data_cpp = Network_Data_Cpp(config_file=sys.argv[1])
network_data_cpp.process_throats()
network_data_cpp.process_pores()
network_data_cpp.process_pore_conns()
#
thrList = network_data_cpp.throat_list
th = network_data_cpp.throat_radius
th = np.array(th) * 2
tl = network_data_cpp.throat_length
tw = network_data_cpp.throat_width
#
conn_in = network_data_cpp.conn_ind_in
conn_out = network_data_cpp.conn_ind_out
#
pc_x = network_data_cpp.pore_coords_x
pc_y = network_data_cpp.pore_coords_y
pc_z = network_data_cpp.pore_coords_z
#
pr = network_data_cpp.pore_radius
#
pl = network_data_cpp.pore_list
p_conn = network_data_cpp.pore_conns
conn_numb = network_data_cpp.conn_number
ppr = network_data_cpp.pore_per_row
#
plx = network_data_cpp.pore_left_x
prx = network_data_cpp.pore_right_x

p_back_y = network_data_cpp.pore_back_y
p_front_y = network_data_cpp.pore_front_y

p_bot_z = network_data_cpp.pore_bot_z
p_top_z = network_data_cpp.pore_top_z

pore_left = p_bot_z
pore_right = p_top_z

langm_coeffs = props.langm_coeff
# # #
# # # # =============================================================================
# # #
hydr_cond = network_data_cpp.hydraulic_cond_coeff

nd_cpp = NetworkDataCpp(thrList, th, tl, tw, conn_in, conn_out, pc_x, pc_y,
                        pc_z, pr, pl, p_conn, conn_numb, ppr, pore_left,
                        pore_right, hydr_cond)
#
eq_pnm = EquationPNM(props_pnm, thrList, th, tl, tw, conn_in, conn_out, pc_x,
                     pc_y, pc_z, pr, pl, p_conn, conn_numb, ppr, pore_left,
                     pore_right, hydr_cond)

diff_pnm = DiffusionPNM(props_pnm, props_diff_vector, thrList, th, tl, tw,
                        conn_in, conn_out, pc_x, pc_y, pc_z, pr, pl, p_conn,
                        conn_numb,
                        ppr, pore_left, pore_right, hydr_cond,
                        langm_coeffs)
# # #
# # # # # =============================================================================
# # # # # Figure 1 (Avg Pore Pressure and Avg Concentration)
# # time = np.arange(0, props.time, props.time_step)
# # # TODO: fix the issue with time step
time = np.linspace(0, props.time, num=int(props.time / props.time_step + 1))
pore_press_av = diff_pnm.get_pressure_av()
matrix_mass_total = diff_pnm.get_matrix_mass_total()
inlet_pressure = diff_pnm.get_inlet_pressure()

matrix_mass_total = matrix_mass_total[:-1]
matrix_mass_total = np.insert(matrix_mass_total, 0, np.nan)
pore_press_av[0] = np.nan
inlet_pressure[0] = np.nan

y_values = {'$M_{matrix}$  ': matrix_mass_total,
            '$P_{av}$': pore_press_av,
            '$P_{in}$': inlet_pressure}

fig_width = 4.5
y_scale = 1.3

# fig, axs = plt.subplots(2, sharex='all',
#                         figsize=(fig_width, fig_width * y_scale),
#                         tight_layout=True)

colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w']

fig, axs = plt.subplots(2, sharex='all', sharey='col',
                        figsize=(fig_width, fig_width * y_scale),
                        tight_layout=True)
plt.suptitle('Model Params vs Time', y=1.0)

# Inflow params (numsticks = 5)
plot_x_ymult(axs[0], time, y_values, 1, 'time, sec', 'mass, $kg$',
             '$P, Pa$', colors, 2, 'solid', [], [])

#
# # Different diffusion
# plot_x_ymult(axs, time, y_values, 1, 'time, sec', 'mass, $kg$',
#              '$P, Pa$', colors, 2, 'solid', [], [])
#
# # df_fig1 = pd.DataFrame({"time": time,
# #                         "Avg_pore_press": pore_press_av,
# #                         "Matrix_vol_total": matrix_mass_total})
# #
# # df_fig1.to_csv(r'../output/fig_press_conc.txt', sep=' ', index=False,
# #                header=True)
# # #
# # # # # =============================================================================
# # # # # Figure 2 (Total Flow Rate)
flow_rate_in = diff_pnm.get_flow_pores_in()
# flow_rate_out = diff_pnm.get_flow_pores_out()
flow_rate_diff = diff_pnm.get_flow_diff()
flow_rate_out = flow_rate_diff + flow_rate_in

flow_rate_out_cum = np.cumsum(flow_rate_out * props.time_step)
flow_rate_diff_cum = np.cumsum(flow_rate_diff * props.time_step)

flow_rate_out_cum = flow_rate_out_cum[:-1]
flow_rate_out_cum = np.insert(flow_rate_out_cum, 0, np.nan)

flow_rate_diff_cum = flow_rate_diff_cum[:-1]
flow_rate_diff_cum = np.insert(flow_rate_diff_cum, 0, np.nan)

flow_rate_diff[0] = np.nan
flow_rate_out[0] = np.nan
#
# # df_fig2 = pd.DataFrame(
# #     {"time": time,
# #      "avg_pore_press": pore_press_av,
# #      "flow_rate_in": flow_rate_in,
# #      "flow_rate_out": flow_rate_out,
# #      "mat_gas_release": flow_rate_diff})
# #
# # df_fig2.to_csv(r'../output/fig_press_flowrates.txt', sep=' ', index=False,
# #                header=True)
# #
y_values = {'$N_{out}$': flow_rate_out_cum,
            '$N_{release}$': flow_rate_diff_cum,
            '$Q_{out}$': flow_rate_out,
            '$Q_{release}$': flow_rate_diff}
#
# # Inflow params (numsticks = 5)
plot_x_ymult(axs[1], time, y_values, 2, 'time, sec', 'mass, $kg$',
             '$Q, kg/sec$', colors, 2, 'solid', [], [])

plt.show()
#
# # Different diffusion
# plot_x_ymult(axs[1], time, y_values, 2, 'time, sec', 'mass, $kg$',
#              '$Q, kg/sec$', colors, 2, 'solid', [], [])
#
# # plt.savefig('../output/flow_params_inflow.eps', format="eps",
# #             bbox_inches='tight')
# plt.show()
#
# # =============================================================================
# # Figure 3 (Langmuir isotherm and density)
# p_av = diff_pnm.get_pressure_av()
#
# # a_dens (kg/m3) is a coefficient for equation density = a * P + b
# a_dens = 6.71079e-06
#
# # b_dens (kg/m3/Pa) is b coefficient for equation density = a * P + b
# b_dens = -2.37253E-02
#
# # rho_const (kg/m3) is a constant density value
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
# y_values = {"Langmuir Conc (kg/m3)": langm_conc_list,
#             "Density (kg/m3)": density_list,
#             "Constant Density (kg/m3)": density_const_list}
#
# plot_x_ymult(p_av, y_values, 1, 'Pressure (Pa)', 'Params vs Pressure')

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

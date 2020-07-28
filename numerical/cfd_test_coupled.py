import sys
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from numerical import PropsPNMCpp
from numerical import EquationPNM
from numerical import Aggregator

from input import Props
from input import Network_Data_Fpnm

from output import plot_x_y
from output import plot_x_ymult
from matplotlib import rc

rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"

# Reading properties in config file
props = Props(config_file=sys.argv[1])
# Getting diffusion properties
props_diff_vector = props.get_diff_props_array()

# Getting pn properties
props_pnm = props.get_props_array()
props_pnm_cpp = PropsPNMCpp(props_pnm)
print(props_pnm_cpp)

# Getting properties of network data
network_data_fpnm = Network_Data_Fpnm(config_file=sys.argv[1])
network_data_fpnm.process_fractures()
network_data_fpnm.process_pores()
#
fractures_list = network_data_fpnm.fractures_list
fractures_widths = network_data_fpnm.fractures_widths
fractures_heights = network_data_fpnm.fractures_heights
fractures_lengths = network_data_fpnm.fractures_lengths

fracs_conn_ind_in = network_data_fpnm.fracs_conn_ind_in
fracs_conn_ind_out = network_data_fpnm.fracs_conn_ind_out

pores_coords_x = network_data_fpnm.pores_coords_x
pores_coords_y = network_data_fpnm.pores_coords_y
pores_coords_z = network_data_fpnm.pores_coords_z

pores_radii = network_data_fpnm.pores_radii
pores_list = network_data_fpnm.pores_list

pores_left_x = network_data_fpnm.pores_left_x
pores_right_x = network_data_fpnm.pores_right_x

pores_back_y = network_data_fpnm.pores_back_y
pores_front_y = network_data_fpnm.pores_front_y

pores_top_z = network_data_fpnm.pores_top_z
pores_bot_z = network_data_fpnm.pores_bot_z

pores_inlet = pores_bot_z
pores_outlet = pores_top_z

langm_coeffs = props.langm_coeff
matrix_volume = props.matrix_volume
solver_method = props.solver_method
# =======================================================================
hydr_cond = network_data_fpnm.hydraulic_cond_coeff
#
eq_pnm = EquationPNM(props_pnm, fractures_list, fractures_heights,
                     fractures_lengths, fractures_widths, fracs_conn_ind_in,
                     fracs_conn_ind_out, pores_coords_x, pores_coords_y,
                     pores_coords_z, pores_radii, pores_list, pores_inlet,
                     pores_outlet, hydr_cond, solver_method)

# eq_pnm.cfd_proc_pure_pnm_dirichlet()

aggregator = Aggregator(props_pnm, props_diff_vector, fractures_list,
                        fractures_heights, fractures_lengths,
                        fractures_widths, fracs_conn_ind_in, fracs_conn_ind_out,
                        pores_coords_x, pores_coords_y, pores_coords_z,
                        pores_radii, pores_list, pores_inlet, pores_outlet,
                        hydr_cond, langm_coeffs, matrix_volume, solver_method)

aggregator.cfd_procedure_pnm_diff()

pore_press_av = aggregator.get_pressure_av()

# =============================================================================
# Figure 1 (Avg Pore Pressure and Avg Concentration)
# time = np.linspace(0, props.time, num=int(props.time / props.time_step + 1))
time = np.cumsum(aggregator.get_time_steps_vec())
time = np.insert(time, 0, 0)
pore_press_av = aggregator.get_pressure_av()
matrix_mass_total = aggregator.get_matrix_mass_total()
inlet_pressure = aggregator.get_inlet_pressure()
#
matrix_mass_total = matrix_mass_total[:-1]
matrix_mass_total = np.insert(matrix_mass_total, 0, np.nan)
pore_press_av[0] = np.nan
inlet_pressure[0] = np.nan
#
y_values = {'$M_{matrix}$': matrix_mass_total,
            '$P_{av}$': pore_press_av,
            '$P_{in}$': inlet_pressure}
#
fig_width = 4.5
y_scale = 1.3
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w']
#
# fig, axs = plt.subplots(2, sharex='all', sharey='col',
#                         figsize=(fig_width, fig_width * y_scale),
#                         tight_layout=True)

fig, axs = plt.subplots(2, sharex='all',
                        figsize=(fig_width, fig_width * y_scale),
                        tight_layout=True)
#
# # Inflow params (numsticks = 5)
plot_x_ymult(axs[0], time, y_values, 1, 'time, sec', 'mass, $kg$',
             '$P, Pa$', colors, 2, 'solid', [6.12 * 10e-8, 6.44 * 10e-8],
             [0.7 + 3e5, 1.7 + 3e5])
#
# # =============================================================================
# # Figure 2 (Total Flow Rate)
#
flow_rate_in = aggregator.get_flow_pores_in()
flow_rate_diff = aggregator.get_flow_diff()
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
y_values = {'$N_{out}$': flow_rate_out_cum,
            '$N_{release}$': flow_rate_diff_cum,
            '$Q_{out}$': flow_rate_out,
            '$Q_{release}$': flow_rate_diff}
# Inflow params (numsticks = 5)
plot_x_ymult(axs[1], time, y_values, 2, 'time, sec', 'mass, $kg$',
             '$Q, kg/sec$', colors, 2, 'solid', [0 * 10e-8, 2.4 * 10e-8],
             [0. * 10e-11, 3.6 * 10e-11])

plt.savefig('../output/flow_params_release.eps', format="eps",
            bbox_inches='tight')
'''
# =============================================================================
# # Figure 3 (Langmuir isotherm and density)
# p_av = diff_pnm.get_pressure_av()
#
# fig_width = 4.5
# y_scale = 0.8
#
# fig, axs = plt.subplots(1, sharex='all', sharey='col',
#                         figsize=(fig_width, fig_width * y_scale),
#                         tight_layout=True)
# # a_dens (kg/m3) is a coefficient for equation density = a * P + b
# a_dens = 6.71079e-06
#
# # b_dens (kg/m3/Pa) is b coefficient for equation density = a * P + b
# b_dens = -2.37253E-02
#
# # rho_const (kg/m3) is a constant density value
# rho_const = a_dens * 300000.3 + b_dens
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
# y_values = {"$C_{langm}$": langm_conc_list,
#             "$\\rho$": density_list,
#             "const $\\rho$": density_const_list}
#
# plot_x_ymult(axs, p_av, y_values, 1, 'Pressure (Pa)',
#              'Concentration, $kg/m^3$', 'Density, $kg/m^3$', colors, 2,
#              'solid', [], [])
#
# plt.show()

# TODO: add one more figure of pressure distribution vs length
# Figure 4 (Pressure and Concentration by Length)
# length_coords = len(pores_coord_x)
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
'''

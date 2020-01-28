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

from output import plot_x_y
from output import plot_x_ymult

props = Props(config_file=sys.argv[1])

props_diff = props.get_diff_props_array()

props_cpp = PropsCpp(props_diff, props.langm_coeff)
print(props_cpp)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++

local_cpp = LocalCpp(props_diff, props.langm_coeff)
convective_cpp = ConvectiveCpp(props_diff, props.langm_coeff)
equation_cpp = EquationCpp(props_diff, props.langm_coeff)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

radius = 4.86487966261438e-06
eff_radius = radius * 5
thr_length = 0.9999723657e-4
conc_wall = 1

local_cpp.calculate_alpha(props_cpp.time, radius, eff_radius, thr_length)

equation_cpp.cfdProcedure(conc_wall, radius, eff_radius, thr_length)

conc = equation_cpp.getConc()
radius_curr = local_cpp.radius_curr
#
flow_rate = equation_cpp.getFlowRate()

grid_centers = []
for i in range(len(radius_curr) - 1):
    grid_centers.append((radius_curr[i + 1] +
                         radius_curr[i]) / 2)

plot_x_y(grid_centers, conc, x_name='Radius (m)',
         y_name='Concentration (kg/m3)',
         graph_name='Concentration distribution,' ' t= 1.1 (sec)')
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print("\n")
#
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

# =============================================================================

hydr_cond = network_data_cpp.hydraulic_cond_coeff

nd_cpp = NetworkDataCpp(thrList, tr, tl, conn_in, conn_out, pc_x, pc_y, pc_z,
                        pr, pl, p_conn, conn_numb, ppr, plx, prx, hydr_cond)

eq_pnm = EquationPNM(props_pnm, thrList, tr, tl, conn_in, conn_out, pc_x, pc_y,
                     pc_z, pr, pl, p_conn, conn_numb, ppr, plx, prx, hydr_cond)

diff_pnm = DiffusionPNM(props_pnm, thrList, tr, tl, conn_in,
                        conn_out, pc_x, pc_y, pc_z, pr, pl, p_conn, conn_numb,
                        ppr, props_diff, props.langm_coeff, plx, prx, hydr_cond)

# =============================================================================
# Figure 1 (Avg Pore Pressure and Avg Concentration)
t = np.arange(0, props.time, props.time_step)
data1 = diff_pnm.get_pressure_av()
data2 = diff_pnm.get_conc_av()

y_values = {"Average Pressure (Pa)": data1, "Avg Concentration (kg/m3)": data2}

plot_x_ymult(t.tolist(), y_values, 1, 'time (sec)', 'FLow Params vs Time')

# =============================================================================
# Figure 2 (Total Flow Rate)
t = np.arange(0, props.time, props.time_step)
data1 = diff_pnm.get_pressure_av()
data2 = diff_pnm.get_flow_pores_out()
data3 = diff_pnm.get_flow_pores_in()

y_values = {"Average Pressure (Pa)": data1, "Outlet Flow Rate (m/sec)": data2,
            "Inlet Flow Rate (m/sec)": data3}

plot_x_ymult(t.tolist(), y_values, 1, 'time (sec)', 'FLow Params vs Time')

# =============================================================================
# Figure 3 (Langmuir isotherm and density)
p_av = diff_pnm.get_pressure_av()

# a_dens (kg/m3) is a coefficient for equation density = a * P + b
a_dens = 6.71079e-06

# b_dens (kg/m3/Pa) is b coefficient for equation density = a * P + b
b_dens = -2.37253E-02

# rho_const (kg/m3) is a constant density value
rho_const = 1.653


def calc_langm_conc(pressure):
    langm_conc = 0
    for i in range(len(props.langm_coeff)):
        langm_conc += props.langm_coeff[i] * pressure ** i
    return langm_conc


langm_conc_list = []
for i in range(len(p_av)):
    langm_conc_list.append(calc_langm_conc(p_av[i]))

density_list = []
for i in range(len(p_av)):
    density_list.append(a_dens * p_av[i] + b_dens)

density_const_list = []
for i in range(len(p_av)):
    density_const_list.append(rho_const)

y_values = {"Langmuir Conc (kg/m3)": langm_conc_list,
            "Density (kg/m3)": density_list,
            "Constant Density (kg/m3)": density_const_list}

plot_x_ymult(p_av, y_values, 1, 'Pressure (Pa)', 'Params vs Pressure')

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

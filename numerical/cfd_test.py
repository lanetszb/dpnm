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

local_cpp = LocalCpp(props_array, props.langm_coeff)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++

convective_cpp = ConvectiveCpp(props_array, props.langm_coeff)

equation_cpp = EquationCpp(props_array, props.langm_coeff)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

radius = 4.86487966261438e-06
eff_radius = radius * 5
thr_length = 0.9999723657e-4
conc_wall = 50

local_cpp.calculate_alpha(props_cpp.time, radius, eff_radius, thr_length)

equation_cpp.cfdProcedure(conc_wall, radius, eff_radius, thr_length)

conc = equation_cpp.getConc()
radius_curr = local_cpp.radius_curr

flow_rate = equation_cpp.getFlowRate()

grid_centers = []
for i in range(len(radius_curr)-1):
    grid_centers.append((radius_curr[i+1] +
                         radius_curr[i])/2)

plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
plt.gca().set(title='Concentration distribution, t = 0.1 (sec)',
              xlabel='Radius (m)',
              ylabel='Concentration (kg/m3)')

print(conc)
plt.plot(grid_centers, conc)

plt.ylim(0, 50)
plt.show()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# print("\n")
# #
# props_pnm = props.get_props_array()
# props_pnm_cpp = PropsPNMCpp(props_pnm)
# print(props_pnm_cpp)
#
# network_data_cpp = Network_Data_Cpp(config_file=sys.argv[1])
# network_data_cpp.process_throats()
# network_data_cpp.process_pores()
# network_data_cpp.process_pore_conns()
# network_data_cpp.process_pore_per_row()
#
# thrList = network_data_cpp.throat_list
# tr = network_data_cpp.throat_radius
# tl = network_data_cpp.throat_length
#
# conn_in = network_data_cpp.conn_ind_in
# conn_out = network_data_cpp.conn_ind_out
#
# pc_x = network_data_cpp.pore_coords_x
# pc_y = network_data_cpp.pore_coords_y
# pc_z = network_data_cpp.pore_coords_z
#
# pr = network_data_cpp.pore_radius
#
# pl = network_data_cpp.pore_list
# p_conn = network_data_cpp.pore_conns
# conn_numb = network_data_cpp.conn_number
# ppr = network_data_cpp.pore_per_row
#
# nd_cpp = NetworkDataCpp(thrList, tr, tl, conn_in, conn_out, pc_x, pc_y, pc_z,
#                         pr, pl, p_conn, conn_numb, ppr)
#
# eq_pnm = EquationPNM(props_pnm, thrList, tr, tl, conn_in, conn_out, pc_x, pc_y,
#                      pc_z, pr, pl, p_conn, conn_numb, ppr)
# #
# diff_pnm = DiffusionPNM(props_pnm, thrList, tr, tl, conn_in,
#                         conn_out, pc_x, pc_y, pc_z, pr, pl, p_conn, conn_numb,
#                         ppr, props_array, props.langm_coeff)

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

from input import Props
from input import Network_Data_Cpp

props = Props(config_file=sys.argv[1])

props_array = props.get_diff_props_array()

props_diff = props.get_props_array()

props_cpp = PropsCpp(props_array)
print(props_cpp)

local_cpp = LocalCpp(props_array)

convective_cpp = ConvectiveCpp(props_array)

equation_cpp = EquationCpp(props_array)

equation_cpp.cfdProcedure(1000)

conc = equation_cpp.getConc()
print(conc)
plt.plot(conc)

plt.ylim(0, 3000)
plt.show()

print("\n")

props_pnm = props.get_props_array()
props_pnm_cpp = PropsPNMCpp(props_pnm)
print(props_pnm_cpp)

network_data_cpp = Network_Data_Cpp(config_file=sys.argv[1])
network_data_cpp.process_throats()
network_data_cpp.process_pores()
network_data_cpp.process_pore_conns()
network_data_cpp.process_pore_per_row()

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

#
nd_cpp = NetworkDataCpp(tr, tl, conn_in, conn_out, pc_x, pc_y, pc_z, pr,
                        pl, p_conn, conn_numb, ppr)

eq_pnm = EquationPNM(props_pnm, tr, tl, conn_in, conn_out, pc_x, pc_y, pc_z, pr,
                     pl, p_conn, conn_numb, ppr)

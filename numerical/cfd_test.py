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

props_pnm_cpp = PropsPNMCpp(props_diff)
print(props_pnm_cpp)

network_data_cpp = Network_Data_Cpp(config_file=sys.argv[1])
network_data_cpp.process_throats()
network_data_cpp.process_pores()

tr = network_data_cpp.pore_radius
tl = network_data_cpp.pore_coords_z

pore_list = network_data_cpp.pore_list

conn_in = network_data_cpp.conn_ind_in
conn_out = network_data_cpp.conn_ind_out

#
nd_cpp = NetworkDataCpp(tr, tl, conn_in, conn_out, tl, conn_in, conn_out, tr,
                        pore_list)

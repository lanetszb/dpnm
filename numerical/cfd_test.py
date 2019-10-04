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

from input import Props

props = Props(config_file=sys.argv[1])
props_array = props.get_diff_props_array()

props_cpp = PropsCpp(props_array)
print(props_cpp)

local_cpp = LocalCpp(props_array)

convective_cpp = ConvectiveCpp(props_array)

equation_cpp = EquationCpp(props_array)

equation_cpp.cfdProcedure(3000)

conc = equation_cpp.getConc()
print(conc)
plt.plot(conc)

plt.ylim(3000, 21000)
plt.show()


import sys
import os
import math

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from numerical import PropsDiffusion
from numerical import EquationDiffusion

from input import Props

from output import plot_x_y
from matplotlib import rc

rc('text', usetex=True)

props = Props(config_file=sys.argv[1])

props_diff_vector = props.get_diff_props_array()
#

# Code for the 1D diffusion only
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Specifying the required parameters
time = props_diff_vector[0]
time_step = props_diff_vector[1]
thr_length = props_diff_vector[2]
radius = props_diff_vector[3]
eff_radius = props_diff_vector[4]
grid_block_n = int(props_diff_vector[5])
conc_fracture_wall = props_diff_vector[6]
diffusivity = props_diff_vector[7]
it_accuracy = props_diff_vector[8]

params_diffusion = {'time': float(time),
                    'time_step': float(time_step),
                    'grid_block_n': int(grid_block_n),
                    'length': float(thr_length),
                    'radius': float(radius),
                    'eff_radius': float(eff_radius),
                    'conc_ini': float(conc_fracture_wall),
                    'diffusivity': float(diffusivity),
                    'it_accuracy': float(it_accuracy)}

props_diff = PropsDiffusion(params_diffusion)
eq_diff = EquationDiffusion(params_diffusion)
eq_diff.cfd_cartesian('dirichlet')

# Evoking parameters required for plotting
conc = eq_diff.conc[1]
L = eff_radius - radius
dX = L / grid_block_n
flow_rate = eq_diff.flow_rate

grid_centers = []
for i in range(grid_block_n):
    grid_centers.append(radius + i * dX + dX / 2)
#
# # Plot conc vs distance
plot_x_y(grid_centers, conc, x_name='Radius (m)',
         y_name='Concentration (kg/m3)',
         graph_name='Concentration distribution',
         line_type='--',
         lw=0.3, marker='s', markersize=4)
# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# # Diffusion Analytical Solution

conc_analyt = []
grid_block_n_analyt = 20
dX = L / grid_block_n_analyt
# #
grid_centers_analyt = []
for i in range(grid_block_n_analyt):
    grid_centers_analyt.append(radius + i * dX + dX / 2)

conc_ini = conc_fracture_wall
conc_out = conc_fracture_wall * 5.0
#
D = diffusivity
t = time
dt = time_step
# #
for i in range(grid_block_n_analyt):
    conc_it = conc_out + (conc_ini - conc_out) * math.erf(
        i * dX / 2 / math.sqrt(D * t))
    conc_analyt.append(conc_it)
# #
conc_analyt.reverse()
#
plot_x_y(grid_centers_analyt, conc_analyt, x_name='Radius (m)',
         y_name='Concentration (kg/m3)',
         graph_name='Concentration distribution',
         line_type='-')

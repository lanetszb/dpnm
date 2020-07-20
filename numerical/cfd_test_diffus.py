import sys
import os
import math

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from numerical import PropsDiffusion
from numerical import LocalDiffusion
from numerical import ConvectiveDiffusion
from numerical import EquationDiffusion

from input import Props

from output import plot_x_y
from matplotlib import rc

rc('text', usetex=True)

props = Props(config_file=sys.argv[1])

props_diff_vector = props.get_diff_props_array()
#
props_diff = PropsDiffusion(props_diff_vector)
print(props_diff, '\n')

# Code for the 1D diffusion only

# Initialising the required classes
convective_diff = ConvectiveDiffusion(props_diff_vector)
local_diff = LocalDiffusion(props_diff_vector)
equation_diff = EquationDiffusion(props_diff_vector)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Specifying the required parameters
time_step = props_diff.time_step
radius = props_diff.radius
eff_radius = props_diff.eff_radius
thr_length = props_diff.length
grid_block_n = props_diff.grid_block_n
conc_fracture_wall = props_diff.conc_ini

# Calculating volumes, both cartesian and cylindrical. Only cylindrical
# calculates radius_curr
# TODO: fix issue with dependence between volumes
# local_diff.calc_vol_cylinder(radius, eff_radius, thr_length)
local_diff.calc_vol_cartesian(radius, eff_radius,
                              thr_length, radius)

# List of FVMs, can be either cylindrical or cartesian
# volume_list = local_diff.vol_cylindr
volume_list = local_diff.vol_cartes

# Calculating omegas, both cartesian and cylindrical.
# TODO: fix issue with dependence between omegas
# convective_diff.calc_omega_cylindr(thr_length, radius, eff_radius)
convective_diff.calc_omega_cartes(radius, thr_length)

# List of omegas, can be either cylindrical or cartesian
# omega_list = convective_diff.omega_cylindr
omega_list = convective_diff.omega_cartes

# Runs the solver, boundary cond can be either:
# 0 (no-flow outer matrix, const C inner) or 1 (Dirichlet, const C both sides).
equation_diff.cfd_procedure('dirichlet', volume_list, omega_list)
# Runs one time step solver used for coupling, boundary conds are mixed
# equation_diff.cfd_procedure_one_step(conc_wall, radius,
#                                      eff_radius, thr_length,
#                                      volume_list, omega_list)

# ==Useless thing used for validation with analytics==
# x = (eff_radius - radius) / (10 - 1)

# equation_diff.cfd_procedure(1, conc_wall, radius - x / 2,
#                             eff_radius + x / 2, thr_length,
#                             volume_list, omega_list)
# ==Useless thing used for validation with analytics==

# Evoking parameters required for plotting
conc = equation_diff.getConc()
local_diff.calc_matr_coord_curr(radius, eff_radius)
radius_curr = local_diff.radius_curr
flow_rate = equation_diff.getFlowRate()

grid_centers = []
for i in range(len(radius_curr) - 1):
    grid_centers.append((radius_curr[i + 1] +
                         radius_curr[i]) / 2)

# Plot conc vs distance
plot_x_y(grid_centers, conc, x_name='Radius (m)',
         y_name='Concentration (kg/m3)',
         graph_name='Concentration distribution',
         line_type='--',
         lw=0.3, marker='s', markersize=4)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Diffusion Analytical Solution
conc_analyt = []
#
# ==Useless thing used for validation with analytics==
# eff_radius_analyt = 2.43243983130719e-05
# radius_analyt = 4.86487966261438e-06
# L = eff_radius_analyt - radius_analyt
# ==Useless thing used for validation with analytics==
#
L = eff_radius - radius

#
grid_block_n_analyt = 299
dX = L / grid_block_n_analyt
#
grid_centers_analyt = []
for i in range(grid_block_n_analyt + 1):
    grid_centers_analyt.append(radius + i * dX)
    # grid_centers_analyt.append(radius + i * dX + dX / 2)
#
conc_ini = conc_fracture_wall
conc_out = conc_fracture_wall * 5.0

D = props_diff.diffusivity
t = props_diff.time
dt = props_diff.time_step
#
for i in range(grid_block_n_analyt + 1):
    conc_it = conc_out + (conc_ini - conc_out) * math.erf(
        i * dX / 2 / math.sqrt(D * t))
    conc_analyt.append(conc_it)
#
conc_analyt.reverse()
#
plot_x_y(grid_centers_analyt, conc_analyt, x_name='Radius (m)',
         y_name='Concentration (kg/m3)',
         graph_name='Concentration distribution',
         line_type='-')

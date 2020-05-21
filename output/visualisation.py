import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from output import plot_x_y
from output import plot_x_ymult
from matplotlib import rc

rc('text', usetex=True)

# Figure validating the diffusion numerical and analytical solutions
table_num = pd.read_csv('conc_time_num.csv')
df_num = pd.DataFrame(table_num)

table_analyt = pd.read_csv('conc_time_analyt.csv')
df_analyt = pd.DataFrame(table_analyt)

fig_width = 4.5
y_scale = 0.9
#
fig1, ax1 = plt.subplots(figsize=(fig_width, fig_width * y_scale),
                         tight_layout=True)

num_arr = [np.array(df_num.coord_num), np.array(df_num.conc_num1),
           np.array(df_num.conc_num2), np.array(df_num.conc_num3)]

for i in range(len(num_arr) - 1):
    plot_x_y(num_arr[0], num_arr[i + 1], x_name='Length, m',
             y_name='Concentration, kg/m$^3$',
             graph_name=None,
             line_type='--',
             lw=1.4,
             color='k')

analyt_arr = [np.array(df_analyt.coord_analyt),
              np.array(df_analyt.conc_analyt1),
              np.array(df_analyt.conc_analyt2),
              np.array(df_analyt.conc_analyt3)]

colors = ['b', 'r', 'g']

for i in range(len(analyt_arr) - 1):
    plot_x_y(analyt_arr[0], analyt_arr[i + 1], x_name='Length, m',
             y_name='Concentration, kg/m$^3$',
             graph_name=None,
             line_type='-',
             color=colors[i],
             alpha=0.5)

line_type = ['--', '-']
line_name = ['numerical', 'analytical']
for i in range(len(line_type)):
    plt.plot([], [], linestyle=line_type[i], c='k',
             label=line_name[i])

legend_1 = plt.legend(scatterpoints=1, frameon=True, labelspacing=1, loc=2)

plt.gca().add_artist(legend_1)

ax1 = plt.scatter([], [], marker="_", c='b', alpha=0.5, s=30,
                  label=str(1))
ax2 = plt.scatter([], [], marker="_", c='r', alpha=0.5, s=30,
                  label=str(7))
ax3 = plt.scatter([], [], marker="_", c='g', alpha=0.5, s=30,
                  label=str(12))

plt.legend(handles=[ax1, ax2, ax3], scatterpoints=1, frameon=True,
           labelspacing=1,
           title='time, sec',
           loc=9)

plt.savefig('../output/diffusion_valid_2.eps', format="eps",
            bbox_inches='tight')

# Figure showing the influence of grid block N on numerical solution

table_gridnum10 = pd.read_csv('gridnum10.csv')
df_gridnum10 = pd.DataFrame(table_gridnum10)

table_gridnum30 = pd.read_csv('gridnum30.csv')
df_gridnum30 = pd.DataFrame(table_gridnum30)

table_gridnum100 = pd.read_csv('gridnum100.csv')
df_gridnum100 = pd.DataFrame(table_gridnum100)

table_gridnum_analyt = pd.read_csv('gridnum_analyt.csv')
df_gridnum_analyt = pd.DataFrame(table_gridnum_analyt)

fig_width = 4.5
y_scale = 0.9
#
fig2, ax2 = plt.subplots(figsize=(fig_width, fig_width * y_scale),
                         tight_layout=True)

coord_analyt = np.array(df_gridnum_analyt.coord_analyt)
conc_analyt = np.array(df_gridnum_analyt.conc_analyt)

plot_x_y(analyt_arr[0], analyt_arr[i + 1], x_name='Length, m',
         y_name='Concentration, kg/m$^3$',
         graph_name=None,
         line_type='-',
         color='k',
         alpha=0.85)

colors = ['g', 'b', 'r']
coord_array_num = [np.array(df_gridnum10.coord_num_10),
                   np.array(df_gridnum30.coord_num_30),
                   np.array(df_gridnum100.coord_num_100)]

conc_array_num = [np.array(df_gridnum10.conc_num_10),
                  np.array(df_gridnum30.conc_num_30),
                  np.array(df_gridnum100.conc_num_100)]

for i in range(len(coord_array_num)):
    plot_x_y(coord_array_num[i], conc_array_num[i], x_name='Length, m',
             y_name='Concentration, kg/m$^3$',
             graph_name=None,
             line_type='--',
             lw=1.4,
             color=colors[i])

line_type = ['--', '-']
line_name = ['numerical', 'analytical']
for i in range(len(line_type)):
    plt.plot([], [], linestyle=line_type[i], c='k',
             label=line_name[i])

legend_1 = plt.legend(scatterpoints=1, frameon=True, labelspacing=1, loc=2)

plt.gca().add_artist(legend_1)

ax1 = plt.scatter([], [], marker="_", c='g', alpha=0.5, s=30,
                  label=str(10))
ax2 = plt.scatter([], [], marker="_", c='b', alpha=0.5, s=30,
                  label=str(30))
ax3 = plt.scatter([], [], marker="_", c='r', alpha=0.5, s=30,
                  label=str(100))

plt.legend(handles=[ax1, ax2, ax3], scatterpoints=1, frameon=True,
           labelspacing=1,
           title='grid number',
           loc=9)

plt.savefig('../output/gridnum_valid.eps', format="eps",
            bbox_inches='tight')
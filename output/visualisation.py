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
import matplotlib.lines as mlines

rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"

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

analyt_arr = [np.array(df_analyt.coord_analyt),
              np.array(df_analyt.conc_analyt1),
              np.array(df_analyt.conc_analyt2),
              np.array(df_analyt.conc_analyt3)]

colors = ['b', 'r', 'g']

for i in range(len(analyt_arr) - 1):
    plot_x_y(analyt_arr[0], analyt_arr[i + 1], x_name='Length, m',
             y_name='Concentration, $kg/m^3$',
             graph_name=None,
             line_type='-',
             color=colors[i],
             lw=1.5)

for i in range(len(num_arr) - 1):
    plot_x_y(num_arr[0], num_arr[i + 1], x_name='Length, $m$',
             y_name='Concentration, $kg/m^3$',
             graph_name=None,
             line_type='-',
             lw=0.3,
             color='k',
             marker='o',
             markersize=2.5)

markers = ['o', '']
lws = [0.3, 1.5]
line_names = ['numerical', 'analytical']
line_types = ['', '-']
for i in range(len(line_names)):
    plt.plot([], [], linestyle=line_types[i], c='k',
             label=line_names[i], marker=markers[i], lw=lws[i],
             markersize='2.5')

legend_1 = plt.legend(scatterpoints=1, frameon=True, labelspacing=1, loc=2)

plt.gca().add_artist(legend_1)

ax1 = mlines.Line2D([], [], linestyle='-', c='b', markersize='2.0',
                    label=str(1))
ax2 = mlines.Line2D([], [], linestyle='-', c='r', markersize='2.0',
                    label=str(7))
ax3 = mlines.Line2D([], [], linestyle='-', c='g', markersize='2.0',
                    label=str(17))

plt.legend(handles=[ax1, ax2, ax3], scatterpoints=1, frameon=True,
           labelspacing=0.5,
           title='time, sec',
           loc=9)

plt.savefig('../output/diffusion_valid_2.eps', format="eps",
            bbox_inches='tight')

# Figure showing the influence of grid block N on numerical solution
table_gridnum10 = pd.read_csv('gridnum10.csv')
df_gridnum10 = pd.DataFrame(table_gridnum10)
#
table_gridnum30 = pd.read_csv('gridnum20.csv')
df_gridnum30 = pd.DataFrame(table_gridnum30)
#
table_gridnum100 = pd.read_csv('gridnum100.csv')
df_gridnum100 = pd.DataFrame(table_gridnum100)
#
table_gridnum_analyt = pd.read_csv('gridnum_analyt.csv')
df_gridnum_analyt = pd.DataFrame(table_gridnum_analyt)
#
fig_width = 4.5
y_scale = 0.9
#
fig2, ax2 = plt.subplots(figsize=(fig_width, fig_width * y_scale),
                         tight_layout=True)
#
coord_analyt = np.array(df_gridnum_analyt.coord_analyt)
conc_analyt = np.array(df_gridnum_analyt.conc_analyt)
#
plot_x_y(analyt_arr[0], analyt_arr[i + 1], x_name='Length, $m$',
         y_name='Concentration, $kg/m^3$',
         graph_name=None,
         line_type='-',
         color='darkorange',
         lw=1.5)

colors = ['g', 'b', 'purple']
#
coord_array_num = [np.array(df_gridnum10.coord_num_10),
                   np.array(df_gridnum30.coord_num_20)]
#
conc_array_num = [np.array(df_gridnum10.conc_num_10),
                  np.array(df_gridnum30.conc_num_20)]
#
markers = ['o', 's']

for i in range(len(coord_array_num)):
    plot_x_y(coord_array_num[i], conc_array_num[i], x_name='Length, $m$',
             y_name='Concentration, $kg/m^3$',
             graph_name=None,
             line_type='-',
             lw=0.7,
             marker=markers[i],
             markersize=3.5,
             color=colors[i])
#
ax1 = mlines.Line2D([], [], linestyle='-', lw=0.7, c='g', marker='o',
                    markersize='3.5',
                    label='10')
ax2 = mlines.Line2D([], [], linestyle='-', lw=0.7, c='b', marker='s',
                    markersize='3.5',
                    label='20')
ax3 = mlines.Line2D([], [], linestyle='-', c='darkorange',
                    label='analytical')
#
legend_2 = plt.legend(handles=[ax3], loc=2)
#
plt.legend(handles=[ax1, ax2], scatterpoints=1,
           frameon=True,
           labelspacing=1,
           title='grid number',
           loc=9)
#
plt.gca().add_artist(legend_2)

plt.savefig('../output/gridnum_valid.eps', format="eps", bbox_inches='tight')

#
# # Figure varying diffusion coefficients

table_vary_diffusion = pd.read_csv('release_vary_diff.csv')
df_vary_diff = pd.DataFrame(table_vary_diffusion)
fig3, ax3 = plt.subplots(figsize=(fig_width, fig_width * y_scale),
                         tight_layout=True)
x_array = np.array(df_vary_diff.time)
diff1_cum = np.array(df_vary_diff.q1cum)
diff2_cum = np.array(df_vary_diff.q2cum)
diff3_cum = np.array(df_vary_diff.q3cum)
diff1 = df_vary_diff.q1
diff2 = df_vary_diff.q2
diff3 = df_vary_diff.q3

y_values = {'$diff1_cum$': diff1_cum, '$diff2_cum$': diff2_cum,
            '$diff3_cum$': diff3_cum,
            '$diff1$': diff1, '$diff2$': diff2,
            '$diff3$': diff3}

colors = ['r', 'g', 'b', 'r', 'g', 'b']
line_type = ['--', '-']
line_name = ['instant', 'accum']
for i in range(len(line_type)):
    plt.plot([], [], linestyle=line_type[i], c='k',
             label=line_name[i])
legend_1 = plt.legend(scatterpoints=1, frameon=True, labelspacing=1,
                      loc='upper center')

plt.gca().add_artist(legend_1)

a1 = plt.scatter([], [], marker="_", c='g', s=30,
                 label=str(5))
a2 = plt.scatter([], [], marker="_", c='b', s=30,
                 label=str(50))
a3 = plt.scatter([], [], marker="_", c='r', s=30,
                 label=str(500))

plt.legend(handles=[a1, a2, a3], scatterpoints=1, frameon=True,
           labelspacing=1, bbox_to_anchor=(0.5, 0.3, 0.5, 0.5),
           title='$D, \; 10^{-11} \: m^2/s$',
           loc='center right')

plot_x_ymult(ax3, x_array, y_values, 3, 'time, sec', 'mass, $kg$',
             '$Q, kg/sec$', colors, 0, 'solid', [0, 2.0 * 10 ** -8],
             [0, 9.5 * 10 ** -11], linestyle='--')

plt.savefig('../output/release_vary_diff.eps', format="eps",
            bbox_inches='tight')

plt.show()

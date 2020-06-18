import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter


def plot_x_y(x_values, y_values, x_name, y_name, graph_name, line_type,
             **kwargs):
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    # plt.suptitle(graph_name, y=1.0, fontsize=17)
    #
    # plt.xlabel(x_name, fontsize=16)
    # plt.ylabel(y_name, fontsize=16)

    plt.suptitle(graph_name)

    plt.xlabel(x_name)
    plt.ylabel(y_name)

    # plt.xlim(min(x_values) - 0.1 * min(x_values),
    #          max(x_values) + 0.1 * max(x_values))
    # plt.ylim(min(y_values), max(y_values) + 0.1 * max(y_values))

    # plt.xlim(1.25e-5, 2.5e-5)
    plt.ylim(min(y_values), max(y_values) + 0.05 * max(y_values))

    plt.plot(x_values, y_values, line_type, **kwargs)

    # plt.show()


def plot_x_ymult(ax, x_values, y_values, y_primary_len, x_name, y1_name,
                 y2_name, colors, marker_size, line_style, y1_lim=[],
                 y2_lim=[], **kwargs):
    colors = colors
    color = 'black'

    ax.set_xlabel(x_name)
    # ax.set_clip_on(False)

    # ax.yaxis.set_major_locator(plt.LinearLocator(numticks=6))
    # ax.yaxis.set_minor_locator(plt.LinearLocator(numticks=5))

    ax.set_ylabel(list(y_values.keys())[0], color=color)
    for idx, value in enumerate((list(y_values.keys())[:y_primary_len])):
        label = list(y_values.keys())[idx]
        # ax.plot(x_values, y_values[value], color=colors[idx], label=label)
        ax.plot(x_values, y_values[value], color=colors[idx], label=label,
                marker='o', markersize=marker_size, linestyle=line_style,
                zorder=10, clip_on=False)
        ax.tick_params(axis='y', labelcolor=color)

    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax.set_ylabel(y1_name)
    if y1_lim:
        ax.set_ylim(y1_lim)

    ax2 = ax.twinx()

    # ax2.yaxis.set_major_locator(plt.LinearLocator(numticks=6))
    # ax2.set_clip_on(False)

    color = 'black'
    ax2.set_ylabel(list(y_values.keys())[y_primary_len], color=color)

    for idx, value in enumerate((list(y_values.keys())[y_primary_len:])):
        label = list(y_values.keys())[y_primary_len + idx]
        # ax2.plot(x_values, y_values[value], color=colors[y_primary_len + idx],
        #          label=label)
        ax2.plot(x_values, y_values[value], color=colors[y_primary_len + idx],
                 label=label, marker='o', markersize=marker_size, **kwargs)
        ax2.tick_params(axis='y', labelcolor=color)

    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylabel(y2_name)
    if y2_lim:
        ax2.set_ylim(y2_lim)

    ax.legend(bbox_to_anchor=(0., -0.2), loc='upper left')
    ax2.legend(bbox_to_anchor=(1., -0.2), loc='upper right')

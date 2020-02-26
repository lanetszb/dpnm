import matplotlib.pyplot as plt


def plot_x_y(x_values, y_values, x_name, y_name, graph_name, line_type):
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    plt.suptitle(graph_name, y=1.0, fontsize=17)

    plt.xlabel(x_name, fontsize=16)
    plt.ylabel(y_name, fontsize=16)

    plt.plot(x_values, y_values, line_type, 'r')
    plt.ylim(min(y_values), max(y_values) + 0.1 * max(y_values))
    plt.show()


def plot_x_ymult(x_values, y_values, y_primary_len, x_name, graph_name):
    fig, ax1 = plt.subplots()
    plt.suptitle(graph_name, y=1.0, fontsize=17)

    colors = ['r', 'g', 'b', 'c']

    color = 'tab:red'
    ax1.set_xlabel(x_name, fontsize=16)
    ax1.set_ylabel(list(y_values.keys())[0], fontsize=16, color=color)

    for idx, value in enumerate((list(y_values.keys())[:y_primary_len])):
        label = list(y_values.keys())[idx]
        ax1.plot(x_values, y_values[value], color=colors[idx], label=label)
        ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()

    color = 'tab:blue'
    ax2.set_ylabel(list(y_values.keys())[y_primary_len], color=color,
                   fontsize=16)

    for idx, value in enumerate((list(y_values.keys())[y_primary_len:])):
        label = list(y_values.keys())[y_primary_len + idx]
        ax2.plot(x_values, y_values[value], color=colors[y_primary_len + idx],
                 label=label)
        ax2.tick_params(axis='y', labelcolor=color)

    ax2.tick_params(axis='y', labelcolor=color)

    fig.legend(loc="upper center", fontsize='x-large')
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()

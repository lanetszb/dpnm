import matplotlib.pyplot as plt


def plot_x_y(x_values, y_values, x_name, y_name, graph_name):
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

    plt.suptitle(graph_name, y=1.0, fontsize=17)

    plt.xlabel(x_name, fontsize=16)
    plt.ylabel(y_name, fontsize=16)

    plt.plot(x_values, y_values)
    plt.ylim(min(y_values), max(y_values) + 0.1 * max(y_values))
    plt.show()

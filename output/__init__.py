import sys
import os

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from output.plot_params import plot_x_y
from output.plot_params import plot_x_ymult

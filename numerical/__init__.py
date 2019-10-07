import sys
import os

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from numerical.cfd import PropsCpp
from numerical.cfd import LocalCpp
from numerical.cfd import ConvectiveCpp
from numerical.cfd import EquationCpp


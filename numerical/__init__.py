import sys
import os

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from numerical.dpnm import PropsDiffusion
from numerical.dpnm import EquationDiffusion
from numerical.dpnm import EquationPNM
from numerical.dpnm import Aggregator

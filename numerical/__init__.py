import sys
import os

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from numerical.cfd import PropsDiffusion
from numerical.cfd import LocalDiffusion
from numerical.cfd import ConvectiveDiffusion
from numerical.cfd import EquationDiffusion
from numerical.cfd import PropsPNMCpp
from numerical.cfd import NetworkDataCpp
from numerical.cfd import EquationPNM
from numerical.cfd import Aggregator

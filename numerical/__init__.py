import sys
import os

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path, '../'))

from numerical.matrix_portrait import Matrix_Portrait

from numerical.cfd import PropsCpp
from numerical.cfd import LocalCpp
from numerical.cfd import ConvectiveCpp
r"""
Package signature_gb for computing noncommutative signature 
Gröbner bases

AUTHORS:

- Clemens Hofstadler (2024-10-28): initial version

"""

# ****************************************************************************
#                          Copyright (C) 2024
#      Clemens Hofstadler(clemens.hofstadler@jku.at)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import pathlib
import sys
sys.path.append(str(pathlib.Path(__file__).parent.resolve()))

from .cython_modules.nc_polynomial import NCPoly
from .cython_modules.labelled_module import LabelledModule
from .cython_modules.matrix_gvw import Matrix_GVW
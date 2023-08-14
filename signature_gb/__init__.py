import pathlib
import sys
sys.path.append(str(pathlib.Path(__file__).parent.resolve()))


from .cython_modules.labelled_module import LabelledModule
from .cython_modules.ncpoly import NCPoly
from .cython_modules.labelled_poly import LabelledPoly
from .cython_modules.labelled_module import LabelledModule
from .cython_modules.sig import Sig
from .python_modules.global_data import *


############################################################################
############################################################################
# Info
############################################################################
############################################################################
print("Package signature_gb version 0.1.0 (beta version)")
print("by Clemens Hofstadler, clemens.hofstadler@mathematik.uni-kassel.de\n")

import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from sage.env import sage_include_directories


ext_module = [
    Extension(
        'cython_modules.linear_algebra',
        sources=['cython_modules/linear_algebra.pyx'],
    ),
    Extension(
        'cython_modules.ncmonomial',
        sources=['cython_modules/ncmonomial.pyx'],
    ),
    Extension(
        'cython_modules.ncpoly',
        sources=['cython_modules/ncpoly.pyx'],
    ),
    Extension(
        'cython_modules.sig',
        sources=['cython_modules/sig.pyx'],
    ),
    Extension(
        'cython_modules.sigpoly',
        sources=['cython_modules/sigpoly.pyx'],
    ),
    Extension(
        'cython_modules.labelled_poly',
        sources=['cython_modules/labelled_poly.pyx'],
    ),
    Extension(
        'cython_modules.orderings',
        sources=['cython_modules/orderings.pyx'],
    ),
    Extension(
        'cython_modules.ambiguity',
        sources=['cython_modules/ambiguity.pyx'],
    ),
    Extension(
        'cython_modules.algorithm',
        sources=['cython_modules/algorithm.pyx'],
    ),
    Extension(
        'cython_modules.sig_gb_algorithm',
        sources=['cython_modules/sig_gb_algorithm.pyx'],
    ),
    Extension(
        'cython_modules.gvw',
        sources=['cython_modules/gvw.pyx'],
    ),
    Extension(
        'cython_modules.matrix_gvw',
        sources=['cython_modules/matrix_gvw.pyx'],
    ),
    Extension(
        'cython_modules.labelled_module',
        sources=['cython_modules/labelled_module.pyx'],
    ),
    Extension(
        'cython_modules.normal_form',
        sources=['cython_modules/normal_form.pyx'],
    ),
]

setup(
  name = 'signature_gb',
  ext_modules = cythonize(ext_module, language_level=3),
  include_dirs = sage_include_directories() + ["."],  
)
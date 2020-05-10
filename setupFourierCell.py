from distutils.core import setup
from Cython.Build import cythonize

setup(name='fourierCell',
    ext_modules=cythonize("fourierCell.pyx", annotate=True))

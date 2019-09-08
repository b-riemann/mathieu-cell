from distutils.core import setup
from Cython.Build import cythonize
from numpy import array

setup(name='Hello world app',
    ext_modules=cythonize("fourierCell.pyx", annotate=True))

from distutils.core import setup
from Cython.Build import cythonize

setup(name='microMacroCell',
      ext_modules=cythonize("microMacroCell.pyx", annotate=True))

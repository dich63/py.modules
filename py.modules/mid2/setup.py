from distutils.core import setup
from Cython.Build import cythonize

setup(name='e_circuit',
      ext_modules=cythonize("e_circuit.py"))
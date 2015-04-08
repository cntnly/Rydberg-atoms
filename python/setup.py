from distutils.core import setup, Extension
from Cython.Build import cythonize

setup(
  name = 'radinte_sqrt',
  ext_modules = cythonize("radinte_sqrt.pyx", library_dirs=['C:\MINGW\lib64']),
)

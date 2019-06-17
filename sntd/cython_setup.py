from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

sourcefiles = ["cython_helpers.pyx"]

extensions = [
    Extension("cython_helpers",
              sourcefiles,
              language="c++",
              include_dirs=[np.get_include()])
]

setup(name="cython_helpers",
      include_dirs=[np.get_include()],
      ext_modules=cythonize(extensions, language_level=3))

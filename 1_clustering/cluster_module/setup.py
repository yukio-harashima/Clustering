# setup.py
# This script compiles the Cython code (.pyx) into a C extension module
# that can be imported into Python.

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

# Define the extension module
extensions = [
    Extension(
        "custom_cluster_fast",  # Name of the compiled module
        ["custom_cluster_fast.pyx"],  # The Cython source file
        include_dirs=[numpy.get_include()],  # Necessary for NumPy C-API
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]
    )
]

setup(
    ext_modules=cythonize(
        extensions,
        compiler_directives={'language_level': "3"} # Use Python 3 syntax
    )
)

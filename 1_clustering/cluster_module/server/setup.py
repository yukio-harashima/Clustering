# setup.py
# This script is used to compile the Cython module (.pyx file).
#
# MODIFICATION:
# - Added compiler and linker flags to enable OpenMP for parallelization.
# - The flags are platform-dependent (GCC/Clang for Linux/macOS, MSVC for Windows).

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import sys

# Define compiler and linker flags for OpenMP
# These flags are different for GCC/Clang (Linux, macOS) and MSVC (Windows)
if sys.platform == 'win32':
    # Windows (MSVC compiler)
    compile_args = ['/openmp']
    link_args = []
else:
    # Linux/macOS (GCC/Clang compilers)
    compile_args = ['-fopenmp']
    link_args = ['-fopenmp']

# Define the extension module
extensions = [
    Extension(
        "custom_cluster_fast",
        ["custom_cluster_fast.pyx"],
        include_dirs=[numpy.get_include()],
        extra_compile_args=compile_args,
        extra_link_args=link_args,
    )
]

# Run the setup
setup(
    ext_modules=cythonize(
        extensions,
        compiler_directives={'language_level': "3"}
    )
)


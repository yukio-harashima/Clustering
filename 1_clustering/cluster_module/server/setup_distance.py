# setup_distance.py
# This script compiles the 'calculate_distances_cython.pyx' module.
# MODIFICATION: Added OpenMP support for parallelization.

import sys
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

# OpenMPフラグをOSに応じて設定
extra_compile_args = []
extra_link_args = []

if sys.platform.startswith('linux') or sys.platform == 'darwin':
    extra_compile_args.append('-fopenmp')
    extra_link_args.append('-fopenmp')
elif sys.platform == 'win32':
    extra_compile_args.append('/openmp')

extensions = [
    Extension(
        "calculate_distances_cython", # ★ 変更点: モジュール名を指定
        ["calculate_distances_cython.pyx"], # ★ 変更点: ソースファイルを指定
        include_dirs=[numpy.get_include()],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args
    )
]

setup(
    ext_modules=cythonize(
        extensions,
        compiler_directives={'language_level': "3"}
    )
)

#!/usr/bin/env python3
from setuptools import setup, Extension
from Cython.Build import cythonize
from glob import glob
import numpy

#-------------Transform Module------------------
filelist1 = ["./EPS09s-Glauber-python-interface/Glauber.pyx", 
			 "./EPS09s-interface-cpp/eps09s.cpp"]
extensions = [
	Extension(
		'Glauber',
		filelist1,
		language="c++",
		extra_compile_args=["-std=c++11"],
		libraries=["m"])
]

setup(
    ext_modules=cythonize(extensions),
    include_dirs=[numpy.get_include()]
        )

#! /usr/bin/env python

# System imports
import platform
import setuptools
from distutils.core import *
from distutils import sysconfig

from distutils.command.build import build

class NativeBuild(build):
    sub_commands = [
        ('build_ext', build.has_ext_modules), 
        ('build_py', build.has_pure_modules),
        ('build_clib', build.has_c_libraries), 
        ('build_scripts', build.has_scripts),
    ]

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

if platform.system().lower() == 'windows':
    _complile_args = [] # ["/EHsc"]
else:
    _complile_args = ["-std=c++11", "-Wno-unused-variable"]

# imagetools extension module
_imagetools = Extension("_imagetools",
                   ["imagetools.i", "cpp/imagetools.cpp", "cpp/geom.cpp", "cpp/raster.cpp",
                       "cpp/csv.cpp", "cpp/ml3d.cpp", "cpp/comparator.cpp"],
                   include_dirs = [numpy_include],
                   swig_opts = ["-c++"],
                   extra_compile_args = _complile_args,
                   )

# imagetools setup
setup(  cmdclass    = {'build': NativeBuild},
        name        = "imagetools",
        description = "Native support for various image segmentation operations",
        author      = "Andrei Volkov",
        version     = "1.0.0",
        license     = "License.txt",
        install_requires=['numpy'],
        ext_modules = [_imagetools,],
        py_modules  = ["imagetools"]
        )

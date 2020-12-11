desc = '''
Linear elastic material model.
'''

from setuptools import setup, Extension

import re
import os
import pybind11
import pyxtensor

header = open('include/GMatElastoPlastic/config.h', 'r').read()
major = re.split(r'(.*)(\#define GMATELASTOPLASTIC_VERSION_MAJOR\ )([0-9]+)(.*)', header)[3]
minor = re.split(r'(.*)(\#define GMATELASTOPLASTIC_VERSION_MINOR\ )([0-9]+)(.*)', header)[3]
patch = re.split(r'(.*)(\#define GMATELASTOPLASTIC_VERSION_PATCH\ )([0-9]+)(.*)', header)[3]

__version__ = '.'.join([major, minor, patch])

include_dirs = [
    os.path.abspath('include/'),
    pyxtensor.find_pyxtensor(),
    pyxtensor.find_pybind11(),
    pyxtensor.find_xtensor(),
    pyxtensor.find_xtl()]

build = pyxtensor.BuildExt

xsimd = pyxtensor.find_xsimd()
if xsimd:
    if len(xsimd) > 0:
        include_dirs += [xsimd]
        build.c_opts['unix'] += ['-march=native', '-DXTENSOR_USE_XSIMD']
        build.c_opts['msvc'] += ['/DXTENSOR_USE_XSIMD']

ext_modules = [Extension(
    'GMatElastoPlastic',
    ['python/main.cpp'],
    include_dirs = include_dirs,
    language = 'c++')]

setup(
    name = 'GMatElastoPlastic',
    description = 'Elasto-plastic material model.',
    long_description = desc,
    keywords = 'Material model; FEM; FFT',
    version = __version__,
    license = 'MIT',
    author = 'Tom de Geus',
    author_email = 'tom@geus.me',
    url = 'https://github.com/tdegeus/GMatElastoPlastic',
    ext_modules = ext_modules,
    install_requires = ['pybind11', 'pyxtensor'],
    cmdclass = {'build_ext': build},
    zip_safe = False)

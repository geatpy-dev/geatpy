# -*- coding: utf-8 -*-
#    This is a list of files to install Geatpy
#
#    Geatpy is a free toolbox: you can redistribute it and/or modify as you want.
#
#    Geatpy is distributed in the hope that it will be useful for the genetic 
#    and evolutionary algorithm, you can get the tutorial from http://www.geatpy.com
#
#    If you want to donate to it, please email geatpy@163.com

import os
import sys
import shutil
import platform
import numpy as np
import setuptools
from Cython.Build import cythonize
from setuptools import setup

LONG_DESCRIPTION = "Geatpy--The Genetic and Evolutionary Algorithms Toolbox for Python. http://www.geatpy.com"

kwargs = dict(name = "geatpy",
    version = "2.2.0",
    description = "Geatpy is a high-performance Genetic and Evolutionary Algorithms toolbox for Python.",
    author = "Geatpy Team",
    author_email = "geatpy@163.com",
    url = "http://www.geatpy.com",
    packages=setuptools.find_packages(),
    include_package_data = True,    # Enabled list file: MANIFEST.in
    install_requires=[
        'numpy>=1.16.0',
        'matplotlib>=3.0.0',
        'scipy>=1.0.0',
    ],
    platforms='any',
    classifiers = [
    'Development Status :: 5 - Production/Stable',
    'Topic :: Software Development',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Scientific/Engineering :: Artificial Intelligence',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    ],
    long_description = LONG_DESCRIPTION,
    zip_safe = False,
)

dest_path = 'geatpy/core/'
source_path = '_core/'

def findAndCopy(path):
    for files in os.listdir(path):
        name = os.path.join(path, files)
        if os.path.isfile(name):
            shutil.copy(name, dest_path + files)

# delete the exist files
def findAndDel(path):
    for files in os.listdir(path):
        name = os.path.join(path, files)
        if os.path.isfile(name):
            os.remove(name)

# copy the core files
if os.path.exists(dest_path) == False:
    os.makedirs(dest_path)
core_path = source_path + platform.system() + '/lib' + platform.architecture()[0][:2] + '/v' + sys.version[:3]
findAndDel(dest_path)
findAndCopy(core_path)

try:
    setup(include_dirs=[np.get_include()],ext_modules=cythonize("build/build.pyx"),language="c", **kwargs),
except Exception:
    setup(**kwargs)

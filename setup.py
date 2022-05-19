# -*- coding: utf-8 -*-
#    This is a list of files to install Geatpy
#
#    Geatpy is a free toolbox: you can redistribute it and/or modify as you want.
#
#    Geatpy is distributed in the hope that it will be useful for the genetic
#    and evolutionary algorithm, you can get the tutorial from http://www.geatpy.com
#
#    If you want to donate to it, please email geatpy@163.com

import pathlib
import platform
import shutil
import sys

import setuptools
from setuptools import setup

dest_path = pathlib.Path('geatpy/core/')
source_path = pathlib.Path('_core/')
if dest_path.exists():
    shutil.rmtree(dest_path)
dest_path.mkdir()
(dest_path / '__init__.py').touch()
python_version = 'v{}.{}'.format(sys.version_info.major,
                                 sys.version_info.minor)
core_path = (source_path / '{}'.format(platform.system())
             / 'lib{}'.format(platform.architecture()[0][:2]) / python_version)
for file in core_path.iterdir():
    shutil.copy(file, dest_path / file.name)
setup(
    name="geatpy",
    version="2.7.0",
    description=("Geatpy is a high-performance "
                 "Genetic and Evolutionary Algorithms toolbox for Python."),
    author="Geatpy Team",
    author_email="geatpy@163.com",
    url="http://www.geatpy.com",
    packages=setuptools.find_packages(),
    include_package_data=True,  # Enabled list file: MANIFEST.in
    python_requires='>=3.5, <3.11',
    install_requires=[
        'numpy>=1.17.0',
        'matplotlib>=3.0.0',
    ],
    platforms='any',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Topic :: Software Development',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Programming Language :: Python',
    ],
    long_description=(
        "Geatpy--The Genetic and Evolutionary Algorithms Toolbox for Python. "
        "http://www.geatpy.com"))

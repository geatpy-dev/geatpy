# -*- coding: utf-8 -*-
#    This is a list of files to install Geatpy
#
#    Geatpy is a free toolbox: you can redistribute it and/or modify as you want.
#
#    Geatpy is distributed in the hope that it will be useful for the genetic 
#    and evolutionary algorithm, you can get the tutorial from http://www.geatpy.com
#
#    If you want to donate to it, please email mathscau@163.com

from setuptools import setup, find_packages

LONG_DESCRIPTION = "Geatpy--The Genetic and Evolutionary Algorithms Toolbox for Python. http://www.geatpy.com"

setup(name = "geatpy",
    version = "2.0.0",
    description = "Geatpy is a high-performance GEA toolbox for Python.",
    author = "Geatpy Team",
    author_email = "mathscau@163.com",
    url = "http://www.geatpy.com",
    packages= find_packages(),
    Scripts = ["geatpy/*"],
    include_package_data = True,    # Enabled list file: MANIFEST.in
    install_requires=[
        'numpy>=1.12.1',
        'matplotlib>=3.0.0',
        'scipy>=1.0.0',
    ],
    classifiers = [
    'Topic :: Utilities',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    ],
    long_description = LONG_DESCRIPTION,
    zip_safe = False,
)

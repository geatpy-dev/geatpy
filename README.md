# **Geatpy** 
The Genetic and Evolutionary Algorithm Toolbox for Python

![Travis](https://travis-ci.org/geatpy-dev/geatpy.svg?branch=master)
![Python](https://img.shields.io/badge/python->=3.5-green.svg)
![Pypi](https://img.shields.io/badge/pypi-2.0.0-blue.svg)

## Introduction
* **Website (including documentation)**: http://www.geatpy.com
* **Demo** : https://github.com/geatpy-dev/geatpy/tree/master/geatpy/demo
* **Pypi page** : https://pypi.org/project/geatpy/
* **Contact us**: http://www.geatpy.com/support
* **Bug reports**: https://github.com/geatpy-dev/geatpy/issues
* **Franchised blog**: https://blog.csdn.net/qq_33353186

It provides:

* global optimization capabilities in **Python** using genetic and other evolutionary algorithms to solve problems unsuitable for traditional optimization approaches.

* a great many of **evolutionary operators**, so that you can deal with **single or multi-objective optimization** problems.

## Installation
1.Installing online:

    pip install geatpy

2.From source:

    python setup.py install

or

    pip install <filename>.whl

**Attention**: Geatpy requires numpy>=1.12.1 and matplotlib>=3.0.0, the installation program won't help you install them so that you have to install both of them by yourselves.

## Versions

**Geatpy** must run under **Python**3.5, 3.6 or 3.7 in x32 or x64 systems.

There are different versions for **Windows**, **Linux** and **Mac**, you can download them from http://geatpy.com/

The version of **Geatpy** on github is the latest version suitable for **Python** >= 3.5

You can also **update** Geatpy by executing the command:

    pip install --upgrade geatpy

If something wrong happened, such as decoding error about 'utf8' of pip, run this command instead or execute it as an administrator:

    pip install --user --upgrade geatpy

Quick start
-----------

Here is the UML of Geatpy2.0.

![image](https://github.com/geatpy-dev/geatpy/blob/master/structure.svg)

You can use **Geatpy** mainly in two steps:

1) Write down the aim function and some relevant settings in a derivative class named **MyProblem**, which is inherited from **Problem** class:

    pip install --upgrade geatpy

2) Instantiate **MyProblem** class and a derivative class inherited from **Algorithm** class in a Python script file "main.py" then execute it. **For example**, trying to find the pareto front of **DTLZ1**, do as the following:

    pip install --upgrade geatpy

The result is:

![image](https://github.com/geatpy-dev/geatpy/blob/master/geatpy/testbed/moea_test/moea_test_DTLZ/Pareto%20Front.png)

The number of non-dominated result: 91

GD:  0.00019492736742063313

IGD:  0.02058320808720775

HV:  0.8413590788841248

Space:  0.00045742613969278813

To get more tutorials, please link to http://www.geatpy.com.

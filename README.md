# Geatpy

The Genetic and Evolutionary Algorithm Toolbox for Python

![Travis](https://www.travis-ci.org/geatpy-dev/geatpy.svg?branch=master)
![Python](https://img.shields.io/badge/python->=3.5-green.svg)
![Pypi](https://img.shields.io/badge/pypi-1.0.0-blue.svg)

Introduction
------------

- **Website (including documentation):** https://www.geatpy.com (repairing)
- **Contact us:** https://www.geatpy.com/supports
- **Source:** https://github.com/geatpy-dev/geatpy
- **Bug reports:** https://github.com/geatpy-dev/geatpy/issues
- **Franchised blog** https://blog.csdn.net/qq_33353186

It provides:

- Global optimization capabilities in Python using genetic and evolutionary algorithm to solve problems unsuitable for traditional optimization approaches.

- A great many of genetic and evolutionary operators, so that you can deal with single or multi-objective optimization problems.

It can work faster with numpy+mkl. If you want to speed your projects, please install numpy+mkl.

Versions
--------------

**Geatpy** must run under **Python**3.5, 3.6 or 3.7 in x32 or x64 in Windows systems.

Installation
------------

1.Via pip::

    pip install geatpy

2.From source::

    python setup.py install
    
    or
    
    pip install <filename>.whl
    
Attention: **Geatpy** requires numpy>=1.10.0 and matplotlib>=1.5.1, the installation program won't help you install them so that you have to install both of them by yourselves.

Quick start
-----------

You can use **Geatpy** mainly in two ways:

1. Create a script, write all the codes on it and run. It's the easiest way, but it needs much too codes and is not good for reuse. To get some examples, please link to http://www.geatpy.com/tutorial.

2. Using templets and functional interfaces. For example, we try to find the Pareto front of DTLZ1, do as the following:

2.1) Write DTLZ1 function on a file named "aimfuc.py" as a functional interfaces:

    """aimfuc.py"""
    # DTLZ1
    def aimfuc(Chrom, M = 3): # M is the dimensions of DTLZ1.
        x = Chrom.T # Chrom is a numpy array standing for the chromosomes of a population
	XM = x[M-1:]
	k = x.shape[0] - M + 1
	gx = 100 * (k + np.sum((XM - 0.5) ** 2 - np.cos(20 * np.pi * (XM - 0.5)), 0))
	ObjV = (np.array([[]]).T) * np.zeros((1, Chrom.shape[0])) # define ObjV to recod function values
	ObjV = np.vstack([ObjV, 0.5 * np.cumprod(x[:M-1], 0)[-1] * (1 + gx)])
	for i in range(2, M):
	    ObjV = np.vstack([ObjV, 0.5 * np.cumprod(x[: M-i], 0)[-1] * (1 - x[M-i]) * (1 + gx)])
	ObjV = np.vstack([ObjV, 0.5 * (1 - x[0]) * (1 + gx)])
	return ObjV.T # use '.T' to change ObjV so that each row stands for function values of each individual of the population

2.2) Write the main script using NSGA-II templet of **Geatpy** to solve the problem.

    """main.py"""
    import numpy as np
    import geatpy as ga # import geatpy
    
    AIM_M = __import__('aimfuc') # get the address of objective function
    
    """==================================variables setting================================"""
    ranges = np.vstack([np.zeros((1,6)), np.ones((1,6))]) # define the ranges of variables in DTLZ1
    borders = np.vstack([np.ones((1,6)), np.ones((1,6))]) # define the borders of variables in DTLZ1
    precisions = [4] * 30 # define the precision of variables in DTLZ1
    """=======================use sga2_templet to find the Pareto front==================="""
    [ObjV, NDSet, times] = ga.nsga2_templet(AIM_M, 'aimfuc',None, None, ranges, borders, precisions, maxormin = 1, MAXGEN = 1000, MAXSIZE = 1000, NIND = 50, SUBPOP = 1, GGAP = 1, selectStyle = 'tour', recombinStyle = 'xovdprs', recopt = 0.9, pm = None, drawing = 1)

The partial of the result is:

![image](https://github.com/geatpy-dev/geatpy/blob/master/geatpy/demo/DTLZ_demo3/Pareto%20Front.png)

To get more tutorials, please link to http://www.geatpy.com

There are more demos in Geatpy's source. Including ZDT1/2/3/4/6、 DTLZ1/2/3/4、single-objective examples、discrete problem solving and so forth.


======
Geatpy
======

The Genetic and Evolutionary Algorithm Toolbox for Python

.. image:: https://www.travis-ci.org/geatpy-dev/geatpy.svg?branch=master
    :target: https://www.travis-ci.org/geatpy-dev/geatpy
    :alt: Build Status

.. image:: https://img.shields.io/badge/python->=3.5-green.svg
    :target: https://www.python.org/downloads/release/python/
    :alt: platform

.. image:: https://img.shields.io/badge/pypi-1.1.4-blue.svg
    :target: https://pypi.org/project/geatpy/
    :alt: versions

Introduction
------------

1. **Website (including documentation):** https://www.geatpy.com
2. **Tutorial pdf:** https://github.com/geatpy-dev/geatpy/tree/master/geatpy/doc/tutorial
3. **Demo:**  https://github.com/geatpy-dev/geatpy/tree/master/geatpy/demo
4. **Contact us:** http://www.geatpy.com/support
5. **Source code:** https://github.com/geatpy-dev/geatpy/tree/master/geatpy/source-code
6. **Bug reports:** https://github.com/geatpy-dev/geatpy/issues
7. **Franchised blog** https://blog.csdn.net/qq_33353186

It provides:

1. global optimization capabilities in Python using genetic and evolutionary algorithm to solve problems unsuitable for traditional optimization approaches.

2. a great many of genetic and evolutionary operators, so that you can deal with single or multi-objective optimization problems.

It can work faster with numpy+mkl. If you want to speed your projects, please install numpy+mkl.

Installation
------------

1. From pip::

    pip install geatpy

2. From source::

    python setup.py install

**Attention**: **Geatpy** requires numpy>=1.12.1 and matplotlib>=2.0.0, the installation program will help you install all the requires. But if something wrong happened, you have to install all requires by yourselves.

Versions
--------------

The version of **Geatpy** on github is the latest version suitable for **Python** >= 3.5

You can also update **Geatpy** by executing the command::

    pip install --upgrade geatpy

Quick start
-----------

You can use **Geatpy** to solve single or multi-objective optimization problems mainly in two ways:

1. Create a script, write all the codes on it and run. It's the easiest way, but it needs much too codes and is not good for reuse. 

2. Using templets and functional interfaces. For example, we try to find the Pareto front of DTLZ1, do as the following:

2.1) Write DTLZ1 function on a file named "aimfuc.py" as a functional interfaces:

    """aimfuc.py"""

    # DTLZ1

    def aimfuc(Chrom, LegV): # LegV is a legal-sign column vector for the population

        M = 3 # M is the dimensions of DTLZ1.
        
        x = Chrom.T # Chrom is a numpy array standing for the chromosomes of the population

        XM = x[M-1:]

        k = x.shape[0] - M + 1

        gx = 100 * (k + np.sum((XM - 0.5) ** 2 - np.cos(20 * np.pi * (XM - 0.5)), 0))

        ObjV = (np.array([[]]).T) * np.zeros((1, Chrom.shape[0])) # define ObjV to recod function values

        ObjV = np.vstack([ObjV, 0.5 * np.cumprod(x[:M-1], 0)[-1] * (1 + gx)])

        for i in range(2, M):

          ObjV = np.vstack([ObjV, 0.5 * np.cumprod(x[: M-i], 0)[-1] * (1 - x[M-i]) * (1 + gx)])

        ObjV = np.vstack([ObjV, 0.5 * (1 - x[0]) * (1 + gx)])

        return [ObjV.T, LegV] # use '.T' to change ObjV so that each row stands for function values of each individual of the population

2.2) Write the main script using NSGA-II templet of **Geatpy** to solve the problem.

    """main.py"""

    import numpy as np

    import geatpy as ga # import geatpy

    AIM_M = __import__('aimfuc') # get the address of objective function DTLZ1
    
    AIM_F = 'DTLZ1' # You can set DTL1,2,3 or 4

    """==================================variables setting================================"""

    ranges = np.vstack([np.zeros((1,7)), np.ones((1,7))]) # define the ranges of variables in DTLZ1
    
    borders = np.vstack([np.ones((1,7)), np.ones((1,7))]) # define the borders of variables in DTLZ1
    
    FieldDR = ga.crtfld(ranges, borders) # create the FieldDR

    """=======================use moea_nsga2_templet to find the Pareto front==================="""

    [ObjV, NDSet, NDSetObjV, times] = ga.nsga2_templet(AIM_M, AIM_F, None, None, FieldDR, problem = 'R', maxormin = 1, MAXGEN = 1000, MAXSIZE = 2000, NIND = 50, SUBPOP = 1, GGAP = 1, selectStyle = 'tour', recombinStyle = 'xovdprs', recopt = 0.9, pm = None, distribute = True, drawing = 1)

The partial of the result can be seen in:

https://github.com/geatpy-dev/geatpy/blob/master/geatpy/demo/DTLZ_demo/Pareto%20Front.png

To get more examples and more tutorials, please link to http://www.geatpy.com/tutorials.

There are also some demos in Geatpy's source. Including ZDT1/2/3/4/6、 DTLZ1/2/3/4、single-objective examples、discrete problem solving and so forth.

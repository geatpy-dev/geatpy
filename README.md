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

1. Write down the aim function and some relevant settings in a derivative class named **MyProblem**, which is inherited from **Problem** class.


    """MyProblem.py"""
    import numpy as np
    import geatpy as ea
    class MyProblem(ea.Problem): # Inherited from Problem class.
        def __init__(self, M):
            self.name = 'DTLZ1' # Problem's name.
            self.M = M # Set the number of objects.
            self.maxormins = [1] * self.M # All objects are need to be minimized.
            self.Dim = self.M + 4 # Set the dimension of decision variables.
            self.varTypes = np.array([0] * self.Dim) # Set the types of decision variables. 0 means continuous while 1 means discrete.
            lb = [0] * self.Dim # The lower bound of each decision variable.
            ub = [1] * self.Dim # The upper bound of each decision variable.
            self.ranges = np.array([lb, ub])
            lbin = [1] * self.Dim # Whether the lower boundary is included.
            ubin = [1] * self.Dim # Whether the upper boundary is included.
            self.borders = np.array([lbin, ubin])
        def aimFuc(self, Vars, CV):
            XM = Vars[:,(self.M-1):]
            g = np.array([100 * (self.Dim - self.M + 1 + np.sum(((XM - 0.5)**2 - np.cos(20 * np.pi * (XM - 0.5))), 1))]).T
            ones_metrix = np.ones((Vars.shape[0], 1))
            ObjV = 0.5 * np.fliplr(np.cumprod(np.hstack([ones_metrix, Vars[:,:self.M-1]]), 1)) * np.hstack([ones_metrix, 1 - Vars[:, range(self.M - 2, -1, -1)]]) * np.tile(1 + g, (1, self.M))
            return ObjV, CV
        def calBest(self):
            uniformPoint, ans = ea.crtup(self.M, 10000) # create 10000 uniform points.
            realBestObjV = uniformPoint / 2
            return realBestObjV

2. Instantiate **MyProblem** class and a derivative class inherited from **Algorithm** class in a Python script file "main.py" then execute it. **For example**, trying to find the pareto front of **DTLZ1**, do as the following:


    """main.py"""
    import geatpy as ea # Import geatpy
    from MyProblem import MyProblem # Import MyProblem class
    """=========================Instantiate your problem=========================="""
    M = 3                      # Set the number of objects.
    problem = MyProblem(M)     # Instantiate MyProblem class
    """===============================Population set=============================="""
    Encoding = 'R'             # Encoding type.
    conordis = 0               # 0 means each element of a chromosome is a continuous number.
    NIND = 100                 # Set the number of individuals.
    Field = ea.crtfld(Encoding, conordis, problem.ranges, problem.borders) # Create the field descriptor.
    population = ea.Population(Encoding, conordis, Field, NIND) # Instantiate Population class(Just instantiate, not initialize the population yet.)
    """================================Algorithm set==============================="""
    myAlgorithm = ea.moea_NSGA3_templet(problem, population) # 实例化一个算法模板对象
    myAlgorithm.MAXGEN = 500 # Set the max times of iteration.
    """===============================Start evolution=============================="""
    NDSet = myAlgorithm.run() # Start evolution.
    """=============================Analyze the result============================="""
    PF = problem.calBest() # Get the global pareto front.
    GD = ea.indicator.GD(NDSet.ObjV, PF) # Calculate GD
    IGD = ea.indicator.IGD(NDSet.ObjV, PF) # Calculate IGD
    HV = ea.indicator.HV(NDSet.ObjV, PF) # Calculate HV
    Space = ea.indicator.spacing(NDSet.ObjV) # Calculate Space
    print('The number of non-dominated result: %s'%(NDSet.sizes))
    print('GD: ',GD)
    print('IGD: ',IGD)
    print('HV: ', HV)
    print('Space: ', Space)

The result is:

![image](https://github.com/geatpy-dev/geatpy/blob/master/geatpy/testbed/moea_test/moea_test_DTLZ/Pareto%20Front.png)

The number of non-dominated result: 91

GD:  0.00019492736742063313

IGD:  0.02058320808720775

HV:  0.8413590788841248

Space:  0.00045742613969278813

To get more tutorials, please link to http://www.geatpy.com.

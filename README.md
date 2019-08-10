# **Geatpy2** 
The Genetic and Evolutionary Algorithm Toolbox for Python

![Travis](https://travis-ci.org/geatpy-dev/geatpy.svg?branch=master)
[![Download](https://img.shields.io/pypi/dm/geatpy.svg)](https://pypi.python.org/pypi/geatpy)
![Python](https://img.shields.io/badge/python->=3.5-green.svg)
![Pypi](https://img.shields.io/badge/pypi-2.2.0-blue.svg)

## Introduction
* **Website (including documentation)**: http://www.geatpy.com
* **Demo** : https://github.com/geatpy-dev/geatpy/tree/master/geatpy/demo
* **Pypi page** : https://pypi.org/project/geatpy/
* **Contact us**: http://geatpy.com/index.php/about/
* **Bug reports**: https://github.com/geatpy-dev/geatpy/issues
* **Notice**: http://geatpy.com/index.php/notice/
* **FAQ**: http://geatpy.com/index.php/faq/

Geatpy provides:

* global optimization capabilities in **Python** using genetic and other evolutionary algorithms to solve problems unsuitable for traditional optimization approaches.

* a great many of **evolutionary operators**, so that you can deal with **single or multi-objective optimization** problems.

## New features of Geatpy 2.2.0

* Add **\<PsyPopulation\> class** and some algorithm templets to support **polysomy population**, which allows to have more than one chromosome in each individual.

* Support evoluating with hybrid encodings.

## Installation
1.Installing online:

    pip install geatpy

2.From source:

    python setup.py install

or

    pip install <filename>.whl

**Attention**: Geatpy requires numpy>=1.16.0, matplotlib>=3.0.0 and scipy>=1.0.0, the installation program won't help you install them so that you have to install both of them by yourselves.

## Versions

**Geatpy** must run under **Python**3.5, 3.6 or 3.7 in Windows x32/x64, Linux x64 or Mac OS x64.

There are different versions for **Windows**, **Linux** and **Mac**, you can download them from http://geatpy.com/

The version of **Geatpy** on github is the latest version suitable for **Python** >= 3.5

You can also **update** Geatpy by executing the command:

    pip install --upgrade geatpy

If something wrong happened, such as decoding error about 'utf8' of pip, run this command instead or execute it as an administrator:

    pip install --user --upgrade geatpy

Quick start
-----------

Here is the UML figure of Geatpy2.

![image](https://github.com/geatpy-dev/geatpy/blob/master/structure.png)

For solving a multi-objective optimization problem, you can use **Geatpy** mainly in two steps:

1.Write down the aim function and some relevant settings in a derivative class named **MyProblem**, which is inherited from **Problem** class:

    """MyProblem.py"""
    import numpy as np
    import geatpy as ea
    class MyProblem(ea.Problem): # Inherited from Problem class.
        def __init__(self, M): # M is the number of objects.
            name = 'DTLZ1' # Problem's name.
            maxormins = [1] * M # All objects are need to be minimized.
            Dim = M + 4 # Set the dimension of decision variables.
            varTypes = [0] * Dim # Set the types of decision variables. 0 means continuous while 1 means discrete.
            lb = [0] * Dim # The lower bound of each decision variable.
            ub = [1] * Dim # The upper bound of each decision variable.
            lbin = [1] * Dim # Whether the lower boundary is included.
            ubin = [1] * Dim # Whether the upper boundary is included.
            # Call the superclass's constructor to complete the instantiation
            ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
        def aimFunc(self, pop): # Write the aim function here, pop is an object of Population class.
            Vars = pop.Phen # Get the decision variables
            XM = Vars[:,(self.M-1):]
            g = np.array([100 * (self.Dim - self.M + 1 + np.sum(((XM - 0.5)**2 - np.cos(20 * np.pi * (XM - 0.5))), 1))]).T
            ones_metrix = np.ones((Vars.shape[0], 1))
            pop.ObjV = 0.5 * np.fliplr(np.cumprod(np.hstack([ones_metrix, Vars[:,:self.M-1]]), 1)) * np.hstack([ones_metrix, 1 - Vars[:, range(self.M - 2, -1, -1)]]) * np.tile(1 + g, (1, self.M))
        def calBest(self): # Calculate the global optimal solution here.
            uniformPoint, ans = ea.crtup(self.M, 10000) # create 10000 uniform points.
            realBestObjV = uniformPoint / 2
            return realBestObjV

2.Instantiate **MyProblem** class and a derivative class inherited from **Algorithm** class in a Python script file "main.py" then execute it. **For example**, trying to find the pareto front of **DTLZ1**, do as the following:

    """main.py"""
    import geatpy as ea # Import geatpy
    from MyProblem import MyProblem # Import MyProblem class
    """=========================Instantiate your problem=========================="""
    M = 3                      # Set the number of objects.
    problem = MyProblem(M)     # Instantiate MyProblem class
    """===============================Population set=============================="""
    Encoding = 'RI'            # Encoding type.
    NIND = 100                 # Set the number of individuals.
    Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders) # Create the field descriptor.
    population = ea.Population(Encoding, Field, NIND) # Instantiate Population class(Just instantiate, not initialize the population yet.)
    """================================Algorithm set==============================="""
    myAlgorithm = ea.moea_NSGA3_templet(problem, population) # Instantiate a algorithm class.
    myAlgorithm.MAXGEN = 500 # Set the max times of iteration.
    """===============================Start evolution=============================="""
    NDSet = myAlgorithm.run() # Run the algorithm templet.
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

Run the "main.py" and the result is:

![image](https://github.com/geatpy-dev/geatpy/blob/master/geatpy/testbed/moea_test/moea_test_DTLZ/Pareto%20Front.svg)

The number of non-dominated result: 91

GD:  0.00019492736742063313

IGD:  0.02058320808720775

HV:  0.8413590788841248

Space:  0.00045742613969278813

For solving another problem: **Ackley-30D**, which has only one object and 30 decision variables, what you need to do is almost the same as above.

1.Write the aim function in "MyProblem.py".

    import numpy as np
    import geatpy as ea
    class Ackley(ea.Problem): # Inherited from Problem class.
        def __init__(self, D = 30):
            name = 'Ackley' # Problem's name.
            M = 1 # Set the number of objects.
            maxormins = [1] * M # All objects are need to be minimized.
            Dim = D # Set the dimension of decision variables.
            varTypes = [0] * Dim # Set the types of decision variables. 0 means continuous while 1 means discrete.
            lb = [-32.768] * Dim # The lower bound of each decision variable.
            ub = [32.768] * Dim # The upper bound of each decision variable.
            lbin = [1] * Dim # Whether the lower boundary is included.
            ubin = [1] * Dim # Whether the upper boundary is included.
            # Call the superclass's constructor to complete the instantiation
            ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
        def aimFunc(self, pop): # Write the aim function here, pop is an object of Population class.
            x = pop.Phen # Get the decision variables
            n = self.Dim
            f = np.array([-20 * np.exp(-0.2*np.sqrt(1/n*np.sum(x**2, 1))) - np.exp(1/n * np.sum(np.cos(2 * np.pi * x), 1)) + np.e + 20]).T
            return f, CV
        def calBest(self): # Calculate the global optimal solution here.
            realBestObjV = np.array([[0]])
            return realBestObjV

2.Write "main.py" to execute the algorithm templet to solve the problem.

    import geatpy as ea # import geatpy
    import numpy as np
    from MyProblem import Ackley
    """=========================Instantiate your problem=========================="""
    problem = Ackley(30) # Instantiate MyProblem class.
    """===============================Population set=============================="""
    Encoding = 'RI'            # Encoding type.
    NIND = 20                  # Set the number of individuals.
    Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders) # Create the field descriptor.
    population = ea.Population(Encoding, Field, NIND) # Instantiate Population class(Just instantiate, not initialize the population yet.)
    """================================Algorithm set==============================="""
    myAlgorithm = ea.soea_DE_rand_1_bin_templet(problem, population) # Instantiate a algorithm class.
    myAlgorithm.MAXGEN = 1000 # Set the max times of iteration.
    myAlgorithm.F = 0.5 # Set the F of DE
    myAlgorithm.pc = 0.2 # Set the Cr of DE (Here it is marked as pc)
    myAlgorithm.drawing = 1 # # 1 means draw the figure of the result
    """===============================Start evolution=============================="""
    [population, obj_trace, var_trace] = myAlgorithm.run() # Run the algorithm templet.
    """=============================Analyze the result============================="""
    best_gen = np.argmin(obj_trace[:, 1]) # Get the best generation.
    best_ObjV = np.min(obj_trace[:, 1])
    print('The objective value of the best solution is: %s'%(best_ObjV))
    print('Effective iteration times: %s'%(obj_trace.shape[0]))
    print('The best generation is: %s'%(best_gen + 1))
    print('The number of evolution is: %s'%(myAlgorithm.evalsNum))

The result is:

![image](https://github.com/geatpy-dev/geatpy/blob/master/geatpy/testbed/soea_test/soea_test_Ackley/result1.svg)

The objective value of the best solution is: 5.8686921988737595e-09

Effective iteration times: 1000

The best generation is: 1000

The number of evolution is: 20000

To get more tutorials, please link to http://www.geatpy.com (Website is maintaining, please visit later!).

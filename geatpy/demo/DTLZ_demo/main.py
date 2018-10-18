# -*- coding: utf-8 -*-
"""
execScript: main.py
description:
    This demo shows how to use Geatpy to test DTLZ-1,2,3,4 multi-objective 
    optimizations' functions.
    The definition of DTLZ-1,2,3,4 functions are written on "aimfuc.py".
    Use 'help' go get more tutorials about the algorithm templet 'moea_nsga2_templet' we use there, 
    or read the source codes in github.
    You can set 'drawing = 2' to get a dynamic graph. (But remember to run the command 'matplotlib qt5' in Python console before.)
"""

import numpy as np
import geatpy as ga # import geatpy

AIM_M = __import__('aimfuc') # get the address of objective function
AIM_F = 'DTLZ1' # You can set DTL1,2,3 or 4

"""==================================variables setting================================"""
ranges = np.vstack([np.zeros((1,7)), np.ones((1,7))]) # define the ranges of variables in DTLZ1
borders = np.vstack([np.ones((1,7)), np.ones((1,7))]) # define the borders of variables in DTLZ1
precisions = [4] * 7 # set any precision that is bigger than 0, or 'crtfld' would create a Field for integer-code.
FieldDR = ga.crtfld(ranges, borders) # create the FieldDR
"""=======================use sga2_templet to find the Pareto front==================="""
[ObjV, NDSet, NDSetObjV, times] = ga.moea_nsga2_templet(AIM_M, AIM_F, None, None, FieldDR, problem = 'R', maxormin = 1, MAXGEN = 500, MAXSIZE = 2000, NIND = 50, SUBPOP = 1, GGAP = 1, selectStyle = 'tour', recombinStyle = 'xovdprs', recopt = 0.9, pm = None, distribute = True, drawing = 1)

"""main.py"""
import numpy as np
import geatpy as ga # import geatpy

AIM_M = __import__('aimfuc') # get the address of objective function
AIM_F = 'DTLZ1' # You can set DTL1,2,3 or 4

"""==================================variables setting================================"""
ranges = np.vstack([np.zeros((1,7)), np.ones((1,7))]) # define the ranges of variables in DTLZ1
borders = np.vstack([np.ones((1,7)), np.ones((1,7))]) # define the borders of variables in DTLZ1
FieldDR = ga.crtfld(ranges, borders) # create the FieldDR
"""=======================use sga2_templet to find the Pareto front==================="""
[ObjV, NDSet, NDSetObjV, times] = ga.nsga2_templet(AIM_M, AIM_F, None, None, FieldDR, problem = 'R', maxormin = 1, MAXGEN = 1000, MAXSIZE = 2000, NIND = 50, SUBPOP = 1, GGAP = 1, selectStyle = 'tour', recombinStyle = 'xovdprs', recopt = 0.9, pm = None, distribute = False, drawing = 2)

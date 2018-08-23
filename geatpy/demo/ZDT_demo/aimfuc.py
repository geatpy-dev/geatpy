# -*- coding: utf-8 -*-
"""
多目标函数ZDT1

Created on Sun Jul 29 14:42:13 2018

@author: jazzbin
"""

import numpy as np

# ZDT1
#def aimfuc(Chrom):
#    
#    ObjV1 = Chrom[:, 0]
#    gx = 1 + (9 / 29) * np.sum(Chrom[:, 1:30], 1)
#    hx = 1 - np.sqrt(ObjV1 / gx)
#    ObjV2 = gx * hx
#    
#    return np.array([ObjV1, ObjV2]).T

# ZDT2
#def aimfuc(Chrom):
#    
#    ObjV1 = Chrom[:, 0]
#    gx = 1 + (9 / 29) * np.sum(Chrom[:, 1:30], 1)
#    hx = 1 - (ObjV1 / gx) ** 2
#    ObjV2 = gx * hx
#    
#    return np.array([ObjV1, ObjV2]).T

# ZDT3
#def aimfuc(Chrom):
#    
#    ObjV1 = Chrom[:, 0]
#    gx = 1 + (9 / 29) * np.sum(Chrom[:, 1:30], 1)
#    hx = 1 - np.sqrt(ObjV1 / gx) - (ObjV1 / gx) * np.sin(10*3.1416*ObjV1)
#    ObjV2 = gx * hx
#    
#    return np.array([ObjV1, ObjV2]).T
    
## ZDT4
#def aimfuc(Chrom):
#    
#    ObjV1 = Chrom[:, 0]
#    gx = 91 + np.sum(Chrom[:, 1:10]**2 - 10 * np.cos(4 * 3.1416 * Chrom[:, 1:10]), 1)
#    hx = 1 - np.sqrt(ObjV1 / gx)
#    ObjV2 = gx * hx
#    
#    return np.array([ObjV1, ObjV2]).T

# ZDT5
def aimfuc(Chrom):
    
    ObjV1 = Chrom[:, 0]
    gx = 1 + (9 / 29) * np.sum(Chrom[:, 1:30], 1)
    hx = 1 - (ObjV1 / gx) ** 2
    ObjV2 = gx * hx
    
    return np.array([ObjV1, ObjV2]).T

# ZDT6
#def aimfuc(Chrom):
#    
#    ObjV1 = 1- np.exp(-4 * Chrom[:, 0]) * (np.sin(6 * 3.1416 * Chrom[:, 0])) ** 6
#    gx = 1 + 9 * (np.sum(Chrom[:, 1:10], 1) / 9) ** 0.25
#    hx = 1 - (ObjV1 / gx) ** 2
#    ObjV2 = gx * hx
#    
#    return np.array([ObjV1, ObjV2]).T
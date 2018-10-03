# -*- coding: utf-8 -*-
"""
aimfc.py - 目标函数demo
描述:
    Geatpy的目标函数遵循本案例的定义方法，
    若要改变目标函数的输入参数、输出参数的格式，则需要修改或自定义算法模板
"""

import numpy as np

def aimfuc(Chrom):
    x1 = Chrom[:, [0]]
    x2 = Chrom[:, [1]]
    x3 = Chrom[:, [2]]
    x4 = Chrom[:, [3]]
    x5 = Chrom[:, [4]]
    ObjV1 = -25 * (x1 - 2)**2 - (x2 - 2)**2 - (x3 - 1)**2 - (x4 - 4)**2 - (x5 - 1)**2
    ObjV2 = (x1 - 1)**2 + (x2 - 1)**2 + (x3 - 1)**2 + (x4 - 1)**2 + (x5 - 1)**2
    
    return np.hstack([ObjV1, ObjV2])
    
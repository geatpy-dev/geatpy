# -*- coding: utf-8 -*-
"""
多目标函数示例1

Created on Sun Jul 29 14:42:13 2018

@author: jazzbin, TRsky
"""

import numpy as np

def aimfuc(Chrom):
    x1 = Chrom[:, 0]
    x2 = Chrom[:, 1]
    x3 = Chrom[:, 2]
    x4 = Chrom[:, 3]
    x5 = Chrom[:, 4]
    ObjV1 = -25 * (x1 - 2)**2 - (x2 - 2)**2 - (x3 - 1)**2 - (x4 - 4)**2 - (x5 - 1)**2
    ObjV2 = (x1 - 1)**2 + (x2 - 1)**2 + (x3 - 1)**2 + (x4 - 1)**2 + (x5 - 1)**2
    
    # 约束条件
    ObjV1[np.where(x1 + x2 < 2)] = ObjV1[np.where(x1 + x2 < 2)] * np.random.rand() / 10
    ObjV2[np.where(x1 + x2 < 2)] = ObjV2[np.where(x1 + x2 < 2)] * (2.1 / np.random.rand())
    
    ObjV1[np.where(x1 + x2 > 6)] = ObjV1[np.where(x1 + x2 > 6)] * np.random.rand() / 10
    ObjV2[np.where(x1 + x2 > 6)] = ObjV2[np.where(x1 + x2 > 6)] * (2.1 / np.random.rand())
    
    ObjV1[np.where(x1 - x2 < -2)] = ObjV1[np.where(x1 - x2 < -2)] * np.random.rand() / 10
    ObjV2[np.where(x1 - x2 < -2)] = ObjV2[np.where(x1 - x2 < -2)] * (2.1 / np.random.rand())
    
    ObjV1[np.where(x1 - 3*x2 > 2)] = ObjV1[np.where(x1 - 3*x2 > 2)] * np.random.rand() / 10
    ObjV2[np.where(x1 - 3*x2 > 2)] = ObjV2[np.where(x1 - 3*x2 > 2)] * (2.1 / np.random.rand())
    
    ObjV1[np.where(4 - (x3 - 3)**2 - x4 < 0)] = ObjV1[np.where(4 - (x3 - 3)**2 - x4 < 0)] * np.random.rand() / 10
    ObjV2[np.where(4 - (x3 - 3)**2 - x4 < 0)] = ObjV2[np.where(4 - (x3 - 3)**2 - x4 < 0)] * (2.1 / np.random.rand())
    
    ObjV1[np.where((x5 - 3)**2 + x4 - 4 < 0)] = ObjV1[np.where((x5 - 3)**2 + x4 - 4 < 0)] * np.random.rand() / 10
    ObjV2[np.where((x5 - 3)**2 + x4 - 4 < 0)] = ObjV2[np.where((x5 - 3)**2 + x4 - 4 < 0)] * (2.1 / np.random.rand())
    
    return np.array([ObjV1, ObjV2]).T
    
    
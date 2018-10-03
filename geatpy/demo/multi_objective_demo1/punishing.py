# -*- coding: utf-8 -*-
"""
punishing.py - 罚函数demo
描述:
    Geatpy的罚函数遵循本案例的定义方法，
    若要改变罚函数的输入参数、输出参数的格式，则需要修改或自定义算法模板
"""

import numpy as np

def punishing(Phen, FitnV):
    x1 = Phen[:, 0]
    x2 = Phen[:, 1]
    x3 = Phen[:, 2]
    x4 = Phen[:, 3]
    x5 = Phen[:, 4]
    # 约束条件
    idx1 = np.where(x1 + x2 < 2)[0]
    idx2 = np.where(x1 + x2 > 6)[0]
    idx3 = np.where(x1 - x2 < -2)[0]
    idx4 = np.where(x1 - 3*x2 > 2)[0]
    idx5 = np.where(4 - (x3 - 3)**2 - x4 < 0)[0]
    idx6 = np.where((x5 - 3)**2 + x4 - 4 < 0)[0]
    FitnV[idx1] = 0
    FitnV[idx2] = 0
    FitnV[idx3] = 0
    FitnV[idx4] = 0
    FitnV[idx5] = 0
    FitnV[idx6] = 0
    exIdx = np.unique(np.hstack([idx1, idx2, idx3, idx4, idx5, idx6])) # 得到非可行解的下标
    return [FitnV, exIdx]
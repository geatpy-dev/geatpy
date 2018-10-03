# -*- coding: utf-8 -*-
"""
punishing.py - 罚函数demo
描述:
    Geatpy的罚函数遵循本案例的定义方法，
    若要改变罚函数的输入参数、输出参数的格式，则需要修改或自定义算法模板
"""

import numpy as np

def punishing(Phen, FitnV):
    x1 = Phen[:, [0]]
    x2 = Phen[:, [1]]
    x3 = 1 - x1 - x2
    # 约束条件
    idx1 = np.where(2 * x1 + x2 > 1)[0]
    idx2 = np.where(x1 + 2 * x3 > 2)[0]
    idx3 = np.where(x1 + x2 > 1)[0]
    # 惩罚
    FitnV[idx1] = 0
    FitnV[idx2] = 0
    FitnV[idx3] = 0
    exIdx = np.unique(np.hstack([idx1, idx2, idx3])) # 得到非可行解在种群中的下标
    return [FitnV, exIdx]

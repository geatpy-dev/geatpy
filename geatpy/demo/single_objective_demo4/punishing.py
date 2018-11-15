# -*- coding: utf-8 -*-
"""
punishing.py - 罚函数demo
描述:
    Geatpy的罚函数遵循本案例的定义方法，
    若要改变罚函数的输入参数、输出参数的格式，则需要修改或自定义算法模板
"""

import numpy as np

def punishing(LegV, FitnV):
    FitnV[np.where(LegV == 0)[0]] = 0
    return FitnV

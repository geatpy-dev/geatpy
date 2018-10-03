# -*- coding: utf-8 -*-
"""
aimfc.py - 目标函数demo
描述:
    Geatpy的目标函数遵循本案例的定义方法，
    若要改变目标函数的输入参数、输出参数的格式，则需要修改或自定义算法模板
"""

import numpy as np

def aimfuc(Phen):
    x1 = Phen[:, [0]]
    x2 = Phen[:, [1]]
    f = 21.5 + x1 * np.sin(4 * np.pi * x1) + x2 * np.sin(20 * np.pi * x2)
    return f

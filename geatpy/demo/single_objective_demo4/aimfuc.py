# -*- coding: utf-8 -*-
"""
aimfc.py - 目标函数demo
描述:
    Geatpy的目标函数遵循本案例的定义方法，
    若要改变目标函数的输入参数、输出参数的格式，则需要修改或自定义算法模板
"""

def aimfuc(Phen):
    x1 = Phen[:, [0]]
    x2 = Phen[:, [1]]
    x3 = 1 - x1 - x2 # 将x1 + x2 + x3 = 1的等式约束降维化处理
    f = 4 * x1 + 2 * x2 + x3
    return f

# -*- coding: utf-8 -*-
"""
执行脚本main.py
描述：
    该demo是展示如何计算带约束的离散型变量的单目标优化问题
    本案例通过调用sga_real_templet算法模板来解决该问题
    其中目标函数写在aimfuc.py文件中
"""

import numpy as np
import geatpy as ga

# 获取函数接口地址
AIM_M = __import__('aimfuc')
# 变量设置
ranges = np.vstack([np.zeros((1, 4)), np.ones((1, 4))]) # 生成自变量的范围矩阵
borders = np.vstack([np.ones((1, 4)), np.ones((1, 4))]) # 生成自变量的边界矩阵
FieldDR = ga.crtfld(ranges, borders) # 生成区域描述器
# 调用编程模板
[pop_trace, var_trace, times] = ga.sga_real_templet(AIM_M, 'aimfuc', None, None, FieldDR, problem = 'I', maxormin = -1, MAXGEN = 50, NIND = 10, SUBPOP = 1, GGAP = 0.9, selectStyle = 'sus', recombinStyle = 'xovdp', recopt = 0.9, pm = 0.1)

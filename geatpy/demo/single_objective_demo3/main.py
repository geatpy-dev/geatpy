# -*- coding: utf-8 -*-
"""
执行脚本main.py
描述：
    该demo是展示如何计算带约束的离散型变量的单目标优化问题
    其中目标函数和约束条件写在aimfuc.py文件中
    max f = 18 * x1 + 10 * x2 + 12 * x3 + 8 * x4
    s.t.
    x1,x2,x3,x4 ∈ {0,1}
    12 * x1 + 6 * x2 + 10 * x3 + 4 * x4 <= 20
    x3 + x4 <= 1
    x3 - x1 <= 0
    x4 - x2 <= 0
    本案例调用了“sga_new_real_templet”这个算法模板，其详细用法可利用help命令查看，或是在github下载并查看源码
    调用算法模板时可以设置drawing=2，此时算法模板将在种群进化过程中绘制动画，但注意执行前要在Python控制台执行命令matplotlib qt5。
"""

import numpy as np
import geatpy as ga

# 获取函数接口地址
AIM_M = __import__('aimfuc')
# 变量设置
ranges = np.vstack([np.zeros((1, 4)), np.ones((1, 4))]) # 生成自变量的范围矩阵
borders = np.vstack([np.ones((1, 4)), np.ones((1, 4))]) # 生成自变量的边界矩阵
FieldDR = ga.crtfld(ranges, borders) # 生成区域描述器
# 调用编程模板(设置problem = 'I'处理离散型变量问题，详见该算法模板的源代码)
[pop_trace, var_trace, times] = ga.sga_new_real_templet(AIM_M, 'aimfuc', None, None, FieldDR, problem = 'I', maxormin = -1, MAXGEN = 50, NIND = 10, SUBPOP = 1, GGAP = 0.9, selectStyle = 'sus', recombinStyle = 'xovdp', recopt = 0.9, pm = 0.1, distribute = True, drawing = 1)

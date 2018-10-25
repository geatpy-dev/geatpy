# -*- coding: utf-8 -*-
"""
描述:
    Rastrigrin函数测试，这是一个多元多峰函数，在取值范围内约有10n个局部极小点，
    是典型的非线性多模态函数，函数图像的峰形高低起伏不定，呈现跳跃性，难以优化找到全局最优解，
    其最小值为：f(0,0,...,0) = 0

Created on Tue Oct 23 22:43:02 2018

@author: TRsky
"""

import numpy as np
import geatpy as ga

# 注意：不建议把目标函数放在执行脚本内，建议放在另一个文件中
def aimfuc(x, LegV): # 定义目标函数
    f = np.array([np.sum((x ** 2 - 10 * np.cos(10 * np.pi * x) + 10), 1)]).T # 调整f使之符合Geatpy的目标函数的数据结构
    return [f, LegV] # 注意返回的参数要符合Geatpy数据结构

if __name__ == "__main__":
    AIM_M = __import__('sga_test_Rastrigrin') # 获取函数所在文件的地址
    # 变量设置
    n = 5
    lb = -5.12 * np.ones((1, n))
    ub = 5.12 * np.ones((1, n))
    ranges = np.vstack([lb, ub])  # 生成自变量的范围矩阵
    lbin = np.ones((1, n))
    ubin = np.ones((1, n))
    borders = np.vstack([lbin, ubin]) # 生成自变量的边界矩阵（1表示变量的区间是闭区间）
    precisions = [1] * n # 根据crtfld的函数特性，这里需要设置精度为任意正值，否则在生成区域描述器时会默认为整数编码，并对变量范围作出一定调整
    FieldDR = ga.crtfld(ranges, borders, precisions) # 生成区域描述器
    # 调用编程模板
    [pop_trace, var_trace, times] = ga.sga_new_real_templet(AIM_M, 'aimfuc', None, None, FieldDR, problem = 'R', maxormin = 1, MAXGEN = 1000, NIND = 100, SUBPOP = 1, GGAP = 0.9, selectStyle = 'sus', recombinStyle = 'xovdprs', recopt = None, pm = None, distribute = True, drawing = 1)

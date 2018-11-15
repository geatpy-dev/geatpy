# -*- coding: utf-8 -*-
"""
描述:
    Pathological病态函数测试，取n=2，这是一个二维复杂函数，具有无数个极小值点，
    有强烈的震荡特性，全局最优解为f(0,0)=0

Created on Tue Oct 23 22:43:02 2018

@author: TRsky
"""

import numpy as np
import geatpy as ga

# 注意：不建议把目标函数放在执行脚本内，建议放在另一个文件中
def aimfuc(Phen, LegV): # 定义目标函数
    x1 = Phen[:, [0]]
    x2 = Phen[:, [1]]
    f = 0.5 + ((np.sin(np.sqrt(x1 ** 2 + x2 ** 2))) ** 2 - 0.5) / ((1 + 0.001 * (x1 ** 2 + x2 ** 2)) ** 2)
    return [f, LegV] # 注意返回的参数要符合Geatpy数据结构

if __name__ == "__main__":
    AIM_M = __import__('sga_test_Pathological') # 获取函数所在文件的地址
    # 变量设置
    ranges = np.array([[-100, -100], [100, 100]])  # 生成自变量的范围矩阵
    borders = np.array([[1, 1], [1, 1]])   # 生成自变量的边界矩阵（1表示变量的区间是闭区间）
    precisions = [1, 1] # 根据crtfld的函数特性，这里需要设置精度为任意正值，否则在生成区域描述器时会默认为整数编码，并对变量范围作出一定调整
    FieldDR = ga.crtfld(ranges, borders, precisions) # 生成区域描述器
    # 调用编程模板
    [pop_trace, var_trace, times] = ga.sga_new_real_templet(AIM_M, 'aimfuc', None, None, FieldDR, problem = 'R', maxormin = 1, MAXGEN = 1000, NIND = 100, SUBPOP = 1, GGAP = 0.9, selectStyle = 'sus', recombinStyle = 'xovdprs', recopt = None, pm = None, distribute = True, drawing = 1)

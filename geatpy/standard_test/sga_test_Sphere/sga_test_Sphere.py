# -*- coding: utf-8 -*-
"""
描述:
    Sphere函数（球函数）测试，取n=3。这是一个简单的平方和函数，只有一个极小点，是多元单峰函数

Created on Tue Oct 23 22:43:02 2018

@author: TRsky
"""

import numpy as np
import geatpy as ga

# 注意：不建议把目标函数放在执行脚本内，建议放在另一个文件中
def aimfuc(x, LegV): # 定义目标函数
    f = np.array([np.sum(x ** 2, 1)]).T # 注意计算结果要符合Geatpy对种群目标函数矩阵的数据结构
    return [f, LegV] # 注意返回的参数要符合Geatpy数据结构

if __name__ == "__main__":
    AIM_M = __import__('sga_test_Sphere') # 获取函数所在文件的地址
    # 变量设置
    ranges = np.array([[-5.12, -5.12, -5.12], [5.12, 5.12, 5.12]])  # 生成自变量的范围矩阵
    borders = np.array([[1, 1, 1], [1, 1, 1]])   # 生成自变量的边界矩阵（1表示变量的区间是闭区间）
    precisions = [1, 1, 1] # 根据crtfld的函数特性，这里需要设置精度为任意正值，否则在生成区域描述器时会默认为整数编码，并对变量范围作出一定调整
    FieldDR = ga.crtfld(ranges, borders, precisions) # 生成区域描述器
    # 调用编程模板
    [pop_trace, var_trace, times] = ga.sga_new_real_templet(AIM_M, 'aimfuc', None, None, FieldDR, problem = 'R', maxormin = 1, MAXGEN = 100, NIND = 50, SUBPOP = 1, GGAP = 0.9, selectStyle = 'sus', recombinStyle = 'xovsp', recopt = None, pm = None, distribute = False, drawing = 1)

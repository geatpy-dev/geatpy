# -*- coding: utf-8 -*-
"""main.py"""
import numpy as np
import geatpy as ga

def aimfuc(x): # define the aim function
    x1 = x[:, 0]; x2 = x[:, 1]
    fun1 = x1**4-10*x1**2+x1*x2+x2**4-x1**2*x2**2
    fun2 = x2**4-x1**2*x2**2+x1**4+x1*x2
    return np.vstack([fun1, fun2]).T # 对矩阵进行转置使得目标函数矩阵符合Geatpy数据结构

if __name__ == "__main__":
    AIM_M = __import__('main') # 获取函数接口所在文件的地址
    # 变量设置
    ranges = np.array([[-5, -5], [5, 5]])  # 生成自变量的范围矩阵
    borders = np.array([[1, 1], [1, 1]])   # 生成自变量的边界矩阵（1表示变量的区间是闭区间）
    precisions = [1, 1] # 因为变量范围都是闭区间，而且nsga2_templet编程模板采用的是实数编码，因此precisions不起作用，随便设置均可
    FieldDR = ga.crtfld(ranges, borders, precisions) # 生成区域描述器
    # 调用编程模板
    [ObjV, NDSet, NDSetObjV, times] = ga.nsga2_templet(AIM_M, 'aimfuc', None, None, FieldDR, 'R', maxormin = 1, MAXGEN = 1000, MAXSIZE = 200, NIND = 200, SUBPOP = 1, GGAP = 1, selectStyle = 'tour', recombinStyle = 'xovdp', recopt = 0.9, pm = 0.6, drawing = 1)
    
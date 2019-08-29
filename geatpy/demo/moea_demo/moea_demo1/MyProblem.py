# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

"""
该案例展示了一个离散决策变量的最小化目标的双目标优化问题。
min f1 = -25 * (x1 - 2)**2 - (x2 - 2)**2 - (x3 - 1)**2 - (x4 - 4)**2 - (x5 - 1)**2
min f2 = (x1 - 1)**2 + (x2 - 1)**2 + (x3 - 1)**2 + (x4 - 1)**2 + (x5 - 1)**2
s.t.
x1 + x2 >= 2
x1 + x2 <= 6
x1 - x2 >= -2
x1 - 3*x2 <= 2
4 - (x3 - 3)**2 - x4 >= 0
(x5 - 3)**2 + x4 - 4 >= 0
x1,x2,x3,x4,x5 ∈ {0,1,2,3,4,5,6,7,8,9,10}
"""

class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self, M = 2):
        name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        Dim = 5 # 初始化Dim（决策变量维数）
        maxormins = [1] * M # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        varTypes = [1] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = [0] * Dim # 决策变量下界
        ub = [10] * Dim # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界
        ubin = [1] * Dim # 决策变量上边界
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        x1 = Vars[:, [0]]
        x2 = Vars[:, [1]]
        x3 = Vars[:, [2]]
        x4 = Vars[:, [3]]
        x5 = Vars[:, [4]]
        f1 = -25 * (x1 - 2)**2 - (x2 - 2)**2 - (x3 - 1)**2 - (x4 - 4)**2 - (x5 - 1)**2
        f2 = (x1 - 1)**2 + (x2 - 1)**2 + (x3 - 1)**2 + (x4 - 1)**2 + (x5 - 1)**2
#        # 利用罚函数法处理约束条件
#        idx1 = np.where(x1 + x2 < 2)[0]
#        idx2 = np.where(x1 + x2 > 6)[0]
#        idx3 = np.where(x1 - x2 < -2)[0]
#        idx4 = np.where(x1 - 3*x2 > 2)[0]
#        idx5 = np.where(4 - (x3 - 3)**2 - x4 < 0)[0]
#        idx6 = np.where((x5 - 3)**2 + x4 - 4 < 0)[0]
#        exIdx = np.unique(np.hstack([idx1, idx2, idx3, idx4, idx5, idx6])) # 得到非可行解的下标
#        f1[exIdx] = f1[exIdx] + np.max(f1) - np.min(f1)
#        f2[exIdx] = f2[exIdx] + np.max(f2) - np.min(f2)
        # 利用可行性法则处理约束条件
        pop.CV = np.hstack([2 - x1 - x2,
                            x1 + x2 - 6,
                            -2 - x1 + x2,
                            x1 - 3*x2 - 2,
                            (x3 - 3)**2 + x4 - 4,
                            4 - (x5 - 3)**2 - x4])
        pop.ObjV = np.hstack([f1, f2]) # 把求得的目标函数值赋值给种群pop的ObjV
    
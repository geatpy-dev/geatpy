# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea
"""
该案例展示了一个带等式约束的连续型决策变量最小化目标的单目标优化问题。
该案例存在多个欺骗性很强的局部最优点，若要增强收敛到全局最优的能力，
如果陷入局部最优，需要调节mutbga变异函数的相关参数，如压缩率、变异梯度等。
max f = 4*x1 + 2*x2 + x3
s.t.
2*x1 + x2 - 1 <= 0
x1 + 2*x3 - 2 <= 0
x1 + x2 + x3 - 1 == 0
0 <= x1,x2,x3 <= 1
"""
class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self):
        self.name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        self.M = 1 # 初始化M（目标维数）
        self.maxormins = [-1] # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = 3 # 初始化Dim（决策变量维数）
        self.varTypes = np.array([1] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [0,0,0] # 决策变量下界
        ub = [1,1,1] # 决策变量上界
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
    
    def aimFuc(self, Vars, CV):
        x1 = Vars[:, [0]]
        x2 = Vars[:, [1]]
        x3 = Vars[:, [2]]
        f = 4*x1 + 2*x2 + x3
        # 采用可行性法则处理约束
        CV = np.hstack([2*x1 + x2 - 1,
                        x1 + 2*x3 - 2,
                        np.abs(x1 + x2 + x3 - 1)])
        return f, CV
    
    def calBest(self):
        realBestObjV = np.array([[2.5]])
        return realBestObjV
    
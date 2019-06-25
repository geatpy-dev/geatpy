# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

"""
该案例展示了一个带约束连续决策变量的最小化目标的双目标优化问题。
min f1 = X**2
min f2 = (X - 2)**2
s.t.
X**2 - 2.5 * X + 1.5 >= 0
10 <= Xi <= 10, (i = 1,2,3,...)
"""

class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self, M = 2):
        self.name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        self.M = M # 初始化M（目标维数）
        self.maxormins = [1] * self.M # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = 1 # 初始化Dim（决策变量维数）
        self.varTypes = np.array([0] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [-10] * self.Dim # 决策变量下界
        ub = [10] * self.Dim # 决策变量上界
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
    
    def aimFuc(self, Vars, CV):
        f1 = Vars**2
        f2 = (Vars - 2)**2
#        # 利用罚函数法处理约束条件
#        exIdx = np.where(Vars**2 - 2.5 * Vars + 1.5 < 0)[0] # 获取不满足约束条件的个体在种群中的下标
#        f1[exIdx] = f1[exIdx] + np.max(f1) - np.min(f1)
#        f2[exIdx] = f2[exIdx] + np.max(f2) - np.min(f2)
        # 利用可行性法则处理约束条件
        CV = -Vars**2 + 2.5 * Vars - 1.5
        
        return np.hstack([f1, f2]), CV
    
    def calBest(self):
        realBestObjV = None
        return realBestObjV
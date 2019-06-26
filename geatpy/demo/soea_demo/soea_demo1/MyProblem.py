# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea
"""
该案例展示了一个简单的连续型决策变量最大化目标的单目标优化问题。
max f = x * np.sin(10 * np.pi * x) + 2.0
s.t.
-1 <= x <= 2
"""
class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self):
        self.name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        self.M = 1 # 初始化M（目标维数）
        self.maxormins = [-1] # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = 1 # 初始化Dim（决策变量维数）
        self.varTypes = np.array([0] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [-1] # 决策变量下界
        ub = [2] # 决策变量上界
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
    
    def aimFuc(self, x, CV):
        f = x * np.sin(10 * np.pi * x) + 2.0
        return f, CV
    
    def calBest(self):
        realBestObjV = None
        return realBestObjV
    
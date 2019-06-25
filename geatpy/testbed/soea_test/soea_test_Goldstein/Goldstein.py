# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ga

class Goldstein(ga.Problem): # 继承Problem父类
    def __init__(self, M = 1, D = 2):
        self.name = 'Goldstein' # 初始化name（函数名称，可以随意设置）
        self.M = M # 初始化M（目标维数）
        self.maxormins = [1] * self.M # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = D # 初始化Dim（决策变量维数）
        self.varTypes = np.array([0] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [-2, -100] # 决策变量下界
        ub = [100, 2] # 决策变量上界
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
    
    def aimFuc(self, Vars, CV):
        x = Vars[:, [0]]
        y = Vars[:, [1]]
        f = (1 + (x + y + 1) ** 2 * (19 - 14 * x + 13 * x ** 2 - 14 * y + 6 * x * y + 3 * y ** 2)) * (30 + (2 * x - 3 * y) ** 2 * (18 - 32 * x + 12 * x ** 2 + 48 * y - 36 * x * y + 27 * y **2))
        return f, CV
    
    def calBest(self):
        realBestObjV = np.array([[3]])
        
        return realBestObjV
    
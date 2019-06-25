# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ga

class Pathological(ga.Problem): # 继承Problem父类
    def __init__(self, M = 1, D = 2):
        self.name = 'Pathological' # 初始化name（函数名称，可以随意设置）
        self.M = M # 初始化M（目标维数）
        self.maxormins = [1] * self.M # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = D # 初始化Dim（决策变量维数）
        self.varTypes = np.array([0] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [-100] * self.Dim # 决策变量下界
        ub = [100] * self.Dim # 决策变量上界
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
    
    def aimFuc(self, Vars, CV):
        x1 = Vars[:, [0]]
        x2 = Vars[:, [1]]
        f = 0.5 + ((np.sin(np.sqrt(x1 ** 2 + x2 ** 2))) ** 2 - 0.5) / ((1 + 0.001 * (x1 ** 2 + x2 ** 2)) ** 2)
        return f, CV
    
    def calBest(self):
        realBestObjV = np.array([[0]])
        
        return realBestObjV
    
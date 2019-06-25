# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class OSY(ea.Problem): # 继承Problem父类
    def __init__(self):
        self.name = 'OSY' # 初始化name（函数名称，可以随意设置）
        self.M = 2 # 初始化M（目标维数）
        self.maxormins = [1] * self.M # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = 6 # 初始化Dim（决策变量维数）
        self.varTypes = np.array([0] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [0, 0, 1, 0, 1, 0] # 决策变量下界全为0
        ub = [10, 10, 5, 6, 5, 10] # 决策变量上界全为1
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
    
    def aimFuc(self, Vars, CV):
        x1 = Vars[:, [0]]
        x2 = Vars[:, [1]]
        x3 = Vars[:, [2]]
        x4 = Vars[:, [3]]
        x5 = Vars[:, [4]]
        x6 = Vars[:, [5]]
        f1 = -(25*(x1-2)**2+(x2-2)**2+(x3-1)**2+(x4-4)**2+(x5-1)**2)
        f2 = x1**2+x2**2+x3**2+x4**2+x5**2+x6**2
        
        # 罚函数法处理约束
#        exIdx1 = np.where(x1+x2-2<0)[0]
#        exIdx2 = np.where(6-x1-x2<0)[0]
#        exIdx3 = np.where(2-x2+x1<0)[0]
#        exIdx4 = np.where(2-x1+3*x2<0)[0]
#        exIdx5 = np.where(4-(x3-3)**2-x4<0)[0]
#        exIdx6 = np.where((x5-3)**2+x6-4<0)[0]
#        exIdx = np.unique(np.hstack([exIdx1, exIdx2, exIdx3, exIdx4, exIdx5, exIdx6]))
#        f1[exIdx] = f1[exIdx] + np.max(f1) - np.min(f1)
#        f2[exIdx] = f2[exIdx] + np.max(f2) - np.min(f2)
        
        # 采用可行性法则处理约束
        CV = np.hstack([-(x1+x2-2),
                       -(6-x1-x2),
                       -(2-x2+x1),
                       -(2-x1+3*x2),
                       -(4-(x3-3)**2-x4),
                       -((x5-3)**2+x6-4)])
        
        return np.hstack([f1, f2]), CV
    
    def calBest(self):
        pass
    
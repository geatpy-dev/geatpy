# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class DTLZ7(ea.Problem): # 继承Problem父类
    def __init__(self, M):
        self.name = 'DTLZ7' # 初始化name（函数名称，可以随意设置）
        self.M = M # 初始化M（目标维数）
        self.maxormins = [1] * self.M # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = self.M + 19 # 初始化Dim（决策变量维数）
        self.varTypes = np.array([0] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [0] * self.Dim # 决策变量下界
        ub = [1] * self.Dim # 决策变量上界
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
    
    def aimFuc(self, Vars, CV):
        _Phen = Vars.copy() # 为了后面计算不影响函数外部的值，这里进行矩阵复制
        ObjV = np.zeros((Vars.shape[0], self.M))
        XM = _Phen[:,(self.M-1):]
        g = np.array([1 + 9 * np.mean(XM, 1)]).T
        ObjV_tmp = Vars[:,:self.M-1]
        ObjV[:,:self.M-1] = ObjV_tmp
        ObjV[:, [self.M-1]] = (1 + g) * (self.M - np.array([np.sum(ObjV_tmp / (1 + np.tile(g, (1, self.M - 1))) * (1 + np.sin(3 * np.pi * ObjV_tmp)), 1)]).T)
        
        return ObjV, CV
    
    def calBest(self):
        realBestObjV = None
        return realBestObjV
    
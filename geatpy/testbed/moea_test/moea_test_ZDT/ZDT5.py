# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class ZDT5(ea.Problem): # 继承Problem父类
    def __init__(self):
        self.name = 'ZDT5' # 初始化name（函数名称，可以随意设置）
        self.M = 2 # 初始化M（目标维数）
        self.maxormins = [1] * self.M # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = 11 # 初始化Dim（决策变量维数）
        self.varTypes = np.array([1] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [0] * self.Dim # 决策变量下界
        ub = [30] + [5] * (self.Dim - 1) # 决策变量上界
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
    
    def aimFuc(self, Vars, CV):
        v = np.zeros(Vars.shape)
        v[np.where(Vars < 5)] += 2
        v[np.where(Vars == 5)] = 1
        ObjV1 = 1 + Vars[:, 0]
        g = np.sum(v[:,1:],1)
        h = 1 / ObjV1
        ObjV2 = g * h
        
        return np.array([ObjV1, ObjV2]).T, CV
    
    def calBest(self):
        ObjV1 = np.array(range(1,32))
        ObjV2 = (self.Dim + 39) / 5 / ObjV1
        realBestObjV = np.array([ObjV1, ObjV2]).T
        
        return realBestObjV
    
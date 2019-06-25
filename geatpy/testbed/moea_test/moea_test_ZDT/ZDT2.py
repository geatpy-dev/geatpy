# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class ZDT2(ea.Problem): # 继承Problem父类
    def __init__(self):
        self.name = 'ZDT2' # 初始化name（函数名称，可以随意设置）
        self.M = 2 # 初始化M（目标维数）
        self.maxormins = [1] * self.M # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = 30 # 初始化Dim（决策变量维数）
        self.varTypes = np.array([0] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [0] * self.Dim # 决策变量下界
        ub = [1] * self.Dim # 决策变量上界
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
    
    def aimFuc(self, Vars, CV):
        ObjV1 = Vars[:, 0]
        gx = 1 + 9 * np.sum(Vars[:, 1:30], 1)
        hx = 1 - (ObjV1 / gx) ** 2
        ObjV2 = gx * hx
        
        return np.array([ObjV1, ObjV2]).T, CV
    
    def calBest(self):
        N = 10000 # 生成10000个参考点
        ObjV1 = np.linspace(0, 1, N)
        ObjV2 = 1 - ObjV1 ** 2
        realBestObjV = np.array([ObjV1, ObjV2]).T
        
        return realBestObjV
    
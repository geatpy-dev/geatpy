# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class DTLZ5(ea.Problem): # 继承Problem父类
    def __init__(self, M):
        self.name = 'DTLZ5' # 初始化name（函数名称，可以随意设置）
        self.M = M # 初始化M（目标维数）
        self.maxormins = [1] * self.M # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = self.M + 9 # 初始化Dim（决策变量维数）
        self.varTypes = np.array([0] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [0] * self.Dim # 决策变量下界
        ub = [1] * self.Dim # 决策变量上界
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
    
    def aimFuc(self, Vars, CV):
        _Phen = Vars.copy() # 为了后面计算不影响函数外部的值，这里进行矩阵复制
        XM = _Phen[:,(self.M-1):]
        g = np.array([np.sum((XM - 0.5) ** 2, 1)]).T
        g_rep   = np.tile(g,(1, self.M - 2)) # 复制g，使其维度与后面与之点乘的矩阵一致
        _Phen[:, 1:self.M-1] = (1 + 2 * g_rep * Vars[:, 1:self.M-1]) / (2 + 2 * g_rep)
        ones_metrix = np.ones((g.shape[0], 1))
        ObjV = np.fliplr(np.cumprod(np.hstack([ones_metrix, np.cos(_Phen[:,:self.M-1] * np.pi / 2)]), 1)) * np.hstack([ones_metrix, np.sin(_Phen[:, range(self.M - 2, -1, -1)] * np.pi / 2)]) * np.tile(1 + g, (1, self.M))
        
        return ObjV, CV
    
    def calBest(self):
        N = 10000 # 生成10000个参考点
        P = np.vstack([np.linspace(0,1,N), np.linspace(1,0,N)]).T
        P = P / np.tile(np.sqrt(np.array([np.sum(P**2, 1)]).T), (1, P.shape[1]))
        P = np.hstack([P[:, np.zeros(self.M-2, dtype = np.int)], P])
        realBestObjV = P / np.sqrt(2) ** np.tile(np.hstack([self.M - 2, np.linspace(self.M - 2, 0, self.M - 1)]), (P.shape[0], 1))
        
        return realBestObjV
    
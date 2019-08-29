# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class ZDT5(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'ZDT5' # 初始化name（函数名称，可以随意设置）
        M = 2 # 初始化M（目标维数）
        maxormins = [1] * M # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 11 # 初始化Dim（决策变量维数）
        varTypes = [1] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = [0] * Dim # 决策变量下界
        ub = [30] + [5] * (Dim - 1) # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界
        ubin = [1] * Dim # 决策变量上边界
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)    
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        v = np.zeros(Vars.shape)
        v[np.where(Vars < 5)] += 2
        v[np.where(Vars == 5)] = 1
        ObjV1 = 1 + Vars[:, 0]
        g = np.sum(v[:,1:],1)
        h = 1 / ObjV1
        ObjV2 = g * h
        pop.ObjV = np.array([ObjV1, ObjV2]).T # 把结果赋值给ObjV
    
    def calBest(self): # 计算全局最优解
        ObjV1 = np.array(range(1,32))
        ObjV2 = (self.Dim + 39) / 5 / ObjV1
        globalBestObjV = np.array([ObjV1, ObjV2]).T
        
        return globalBestObjV
    
# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class ZDT6(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'ZDT6' # 初始化name（函数名称，可以随意设置）
        M = 2 # 初始化M（目标维数）
        maxormins = [1] * M # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 10 # 初始化Dim（决策变量维数）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = [0] * Dim # 决策变量下界
        ub = [1] * Dim # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界
        ubin = [1] * Dim # 决策变量上边界
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        ObjV1 = 1 - np.exp(-4 * Vars[:, 0]) * (np.sin(6 * np.pi * Vars[:, 0])) ** 6
        gx = 1 + 9 * (np.sum(Vars[:, 1:10], 1) / 9) ** 0.25
        hx = 1 - (ObjV1 / gx) ** 2
        ObjV2 = gx * hx
        pop.ObjV = np.array([ObjV1, ObjV2]).T # 把结果赋值给ObjV
    
    def calBest(self): # 计算全局最优解
        N = 10000 # 生成10000个参考点
        ObjV1 = np.linspace(0.280775, 1, N)
        ObjV2 = 1 - ObjV1 ** 2;
        globalBestObjV = np.array([ObjV1, ObjV2]).T
        
        return globalBestObjV
    
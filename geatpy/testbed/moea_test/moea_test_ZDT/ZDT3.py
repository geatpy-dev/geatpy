# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class ZDT3(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'ZDT3' # 初始化name（函数名称，可以随意设置）
        M = 2 # 初始化M（目标维数）
        maxormins = [1] * M # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 30 # 初始化Dim（决策变量维数）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = [0] * Dim # 决策变量下界
        ub = [1] * Dim # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界
        ubin = [1] * Dim # 决策变量上边界
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        ObjV1 = Vars[:, 0]
        gx = 1 + 9 * np.sum(Vars[:, 1:30], 1)
        hx = 1 - np.sqrt(ObjV1 / gx) - (ObjV1 / gx) * np.sin(10*3.1416*ObjV1)
        ObjV2 = gx * hx
        pop.ObjV = np.array([ObjV1, ObjV2]).T # 把结果赋值给ObjV
    
    def calBest(self): # 计算全局最优解
        N = 10000 # 生成10000个参考点
        ObjV1 = np.linspace(0, 1, N)
        ObjV2 = 1 - ObjV1**0.5 - ObjV1 * np.sin(10*np.pi * ObjV1)
        f = np.array([ObjV1, ObjV2]).T
        levels, criLevel = ea.ndsortESS(f, None, 1)
        globalBestObjV = f[np.where(levels == 1)[0]]
        
        return globalBestObjV
    
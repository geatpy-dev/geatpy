# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class Goldstein(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'Goldstein' # 初始化name（函数名称，可以随意设置）
        M = 1 # 初始化M（目标维数）
        Dim = 2 # 初始化Dim（决策变量维数）
        maxormins = [1] * M # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = [-2, -100] # 决策变量下界
        ub = [100, 2] # 决策变量上界
        lbin = [1] * Dim
        ubin = [1] * Dim
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        x = Vars[:, [0]]
        y = Vars[:, [1]]
        pop.ObjV = (1 + (x + y + 1) ** 2 * (19 - 14 * x + 13 * x ** 2 - 14 * y + 6 * x * y + 3 * y ** 2)) * (30 + (2 * x - 3 * y) ** 2 * (18 - 32 * x + 12 * x ** 2 + 48 * y - 36 * x * y + 27 * y **2))
    
    def calBest(self): # 计算全局最优解
        globalBestObjV = np.array([[3]])
        return globalBestObjV
    
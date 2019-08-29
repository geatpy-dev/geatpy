# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class Fonseca(ea.Problem): # 继承Problem父类
    def __init__(self, M = 2, Dim = 3):
        name = 'Fonseca' # 初始化name（函数名称，可以随意设置）
        maxormins = [1] * M # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = [-4] * Dim # 决策变量下界
        ub = [4] * Dim # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界
        ubin = [1] * Dim # 决策变量上边界
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        f1 = 1 - np.exp(-np.sum((Vars - 1 / np.sqrt(self.Dim))**2, 1, keepdims = True))
        f2 = 1 - np.exp(-np.sum((Vars + 1 / np.sqrt(self.Dim))**2, 1, keepdims = True))
        pop.ObjV = np.hstack([f1, f2]) # 把求得的目标函数值赋值给种群pop的ObjV
    
    def calBest(self): # 计算全局最优解
        return None
    
# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class Rosenbrock(ea.Problem): # 继承Problem父类
    def __init__(self, Dim = 2):
        name = 'Rosenbrock' # 初始化name（函数名称，可以随意设置）
        M = 1 # 初始化M（目标维数）
        maxormins = [1] * M # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = [-2.048] * Dim # 决策变量下界
        ub = [2.048] * Dim # 决策变量上界
        lbin = [1] * Dim
        ubin = [1] * Dim
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    
    def aimFunc(self, pop): # 目标函数
        x = pop.Phen # 得到决策变量矩阵
        Nvar = self.Dim
        Mat1 = x[:, :Nvar - 1]
        Mat2 = x[:, 1:Nvar]
        pop.ObjV = np.array([np.sum((100*(Mat2 - Mat1**2)**2 + (1 - Mat1)**2).T, 0)]).T

    def calBest(self): # 计算全局最优解
        globalBestObjV = np.array([[0]])
        return globalBestObjV
    
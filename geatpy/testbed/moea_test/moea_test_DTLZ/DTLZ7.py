# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class DTLZ7(ea.Problem): # 继承Problem父类
    def __init__(self, M):
        name = 'DTLZ7' # 初始化name（函数名称，可以随意设置）
        maxormins = [1] * M # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = M + 19 # 初始化Dim（决策变量维数）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = [0] * Dim # 决策变量下界
        ub = [1] * Dim # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界
        ubin = [1] * Dim # 决策变量上边界
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        _Phen = Vars.copy() # 为了后面计算不影响函数外部的值，这里进行矩阵复制
        ObjV = np.zeros((Vars.shape[0], self.M))
        XM = _Phen[:,(self.M-1):]
        g = 1 + 9 * np.mean(XM, 1, keepdims = True)
        ObjV_tmp = Vars[:,:self.M-1]
        ObjV[:,:self.M-1] = ObjV_tmp
        ObjV[:, [self.M-1]] = (1 + g) * (self.M - np.sum(ObjV_tmp / (1 + np.tile(g, (1, self.M - 1))) * (1 + np.sin(3 * np.pi * ObjV_tmp)), 1, keepdims = True))
        pop.ObjV = ObjV # 把求得的目标函数值赋值给种群pop的ObjV
    
    def calBest(self): # 计算全局最优解
        print('暂无DTLZ-7的真实前沿数据。')
        globalBestObjV = None
        return globalBestObjV
    
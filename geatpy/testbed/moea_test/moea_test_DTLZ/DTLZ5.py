# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class DTLZ5(ea.Problem): # 继承Problem父类
    def __init__(self, M):
        name = 'DTLZ5' # 初始化name（函数名称，可以随意设置）
        maxormins = [1] * M # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = M + 9 # 初始化Dim（决策变量维数）
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
        XM = _Phen[:,(self.M-1):]
        g = np.sum((XM - 0.5) ** 2, 1, keepdims = True)
        g_rep   = np.tile(g,(1, self.M - 2)) # 复制g，使其维度与后面与之点乘的矩阵一致
        _Phen[:, 1:self.M-1] = (1 + 2 * g_rep * Vars[:, 1:self.M-1]) / (2 + 2 * g_rep)
        ones_metrix = np.ones((g.shape[0], 1))
        f = np.fliplr(np.cumprod(np.hstack([ones_metrix, np.cos(_Phen[:,:self.M-1] * np.pi / 2)]), 1)) * np.hstack([ones_metrix, np.sin(_Phen[:, range(self.M - 2, -1, -1)] * np.pi / 2)]) * np.tile(1 + g, (1, self.M))
        pop.ObjV = f # 把求得的目标函数值赋值给种群pop的ObjV
    
    def calBest(self): # 计算全局最优解
        N = 10000 # 生成10000个参考点
        P = np.vstack([np.linspace(0,1,N), np.linspace(1,0,N)]).T
        P = P / np.tile(np.sqrt(np.sum(P**2, 1, keepdims = True)), (1, P.shape[1]))
        P = np.hstack([P[:, np.zeros(self.M-2, dtype = np.int)], P])
        globalBestObjV = P / np.sqrt(2) ** np.tile(np.hstack([self.M - 2, np.linspace(self.M - 2, 0, self.M - 1)]), (P.shape[0], 1))
        
        return globalBestObjV
    
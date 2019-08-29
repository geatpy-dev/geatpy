# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class C2_DTLZ2(ea.Problem): # 继承Problem父类
    def __init__(self, M = 3):
        name = 'C2-DTLZ2' # 初始化name（函数名称，可以随意设置）
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
        XM = Vars[:,(self.M-1):]
        g = np.sum((XM - 0.5) ** 2, 1, keepdims = True)
        ones_metrix = np.ones((g.shape[0], 1))
        f = np.fliplr(np.cumprod(np.hstack([ones_metrix, np.cos(Vars[:,:self.M-1] * np.pi / 2)]), 1)) * np.hstack([ones_metrix, np.sin(Vars[:, range(self.M - 2, -1, -1)] * np.pi / 2)]) * np.tile(1 + g, (1, self.M))
        # 计算违反约束程度矩阵的值
        r = 0.4 if self.M == 3 else 0.5
        CV = np.min([np.min((f-1)**2 + np.sum(f**2, 1, keepdims = True) - f**2 - r**2, 1, keepdims = True), np.sum((f-1/np.sqrt(self.M))**2, 1, keepdims = True) - r**2], 0)
        pop.ObjV = f # 把求得的目标函数值赋值给种群pop的ObjV
        pop.CV = CV # 把求得的违反约束程度矩阵赋值给种群pop的CV
    
    def calBest(self): # 计算全局最优解
        Point, ans = ea.crtup(self.M, 10000) # 生成10000个在各目标的单位维度上均匀分布的参考点
        Point /= np.tile(np.sqrt(np.sum(Point**2, 1, keepdims = True)), (1, self.M))
        r = 0.4 if self.M == 3 else 0.5
        globalBestObjV = Point[np.where(np.min([np.min((Point-1)**2 + np.sum(Point**2, 1, keepdims = True) - Point**2-r**2, 1, keepdims = True), np.sum((Point-1/np.sqrt(self.M))**2, 1, keepdims = True) - r**2], 0) <= 0)[0], :]
        return globalBestObjV
    
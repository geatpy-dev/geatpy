# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class OSY(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'OSY' # 初始化name（函数名称，可以随意设置）
        M = 2 # 初始化M（目标维数）
        maxormins = [1] * M # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 6 # 初始化Dim（决策变量维数）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = [0, 0, 1, 0, 1, 0] # 决策变量下界全为0
        ub = [10, 10, 5, 6, 5, 10] # 决策变量上界全为1
        lbin = [1] * Dim # 决策变量下边界
        ubin = [1] * Dim # 决策变量上边界
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
        
    def aimFunc(self, pop): # 目标函数
        X = pop.Phen # 得到决策变量矩阵
        x1 = X[:, [0]]
        x2 = X[:, [1]]
        x3 = X[:, [2]]
        x4 = X[:, [3]]
        x5 = X[:, [4]]
        x6 = X[:, [5]]
        f1 = -(25*(x1-2)**2+(x2-2)**2+(x3-1)**2+(x4-4)**2+(x5-1)**2)
        f2 = np.sum(X**2, 1, keepdims = True)
#        # 罚函数法处理约束
#        exIdx1 = np.where(x1+x2-2<0)[0]
#        exIdx2 = np.where(6-x1-x2<0)[0]
#        exIdx3 = np.where(2-x2+x1<0)[0]
#        exIdx4 = np.where(2-x1+3*x2<0)[0]
#        exIdx5 = np.where(4-(x3-3)**2-x4<0)[0]
#        exIdx6 = np.where((x5-3)**2+x6-4<0)[0]
#        exIdx = np.unique(np.hstack([exIdx1, exIdx2, exIdx3, exIdx4, exIdx5, exIdx6]))
#        alpha = 2 # 惩罚缩放因子
#        beta = 1 # 惩罚最小偏移量
#        f1[exIdx] = f1[exIdx] + self.maxormins[0] * alpha * (np.max(f1)-np.min(f1)+beta)
#        f2[exIdx] = f2[exIdx] + self.maxormins[1] * alpha * (np.max(f2)-np.min(f2)+beta)
        # 采用可行性法则处理约束
        pop.CV = np.hstack([-(x1+x2-2)/2,
                       -(6-x1-x2)/6,
                       -(2-x2+x1)/2,
                       -(2-x1+3*x2)/2,
                       -(4-(x3-3)**2-x4)/4,
                       -((x5-3)**2+x6-4)/4])
        
        pop.ObjV = np.hstack([f1, f2]) # 把求得的目标函数值赋值给种群pop的ObjV
    
    
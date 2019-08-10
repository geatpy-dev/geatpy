# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class WFG3(ea.Problem): # 继承Problem父类
    def __init__(self, M = 3):
        name = 'WFG3' # 初始化name（函数名称，可以随意设置）
        maxormins = [1] * M # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = M + 9 # 初始化Dim（决策变量维数）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = [0] * Dim # 决策变量下界
        ub = list(range(2, 2 * Dim + 1, 2)) # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界
        ubin = [1] * Dim # 决策变量上边界
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
        # 目标函数中用到的一些参数设置
        self.K = self.M - 1
        self.L = self.Dim - self.K
        self.S = np.array([list(range(2, 2 * self.M + 1, 2))])
        self.D = 1
        self.A = [1] + [0] * (self.M - 2)
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        N, Lind = Vars.shape
        M = self.M
        K = self.K
        L = self.L
        S = self.S
        D = self.D
        A = np.tile(self.A, (N, 1))
        Z = Vars / np.tile(np.array([range(2, Lind * 2 + 1, 2)]), (N,1))
        t1 = np.zeros((N, K + L))
        t1[:, :K] = Z[:,:K]
        t1[:, K:] = s_linear(Z[:,K:], 0.35)
        t2 = np.zeros((N, int(K + L / 2)))
        t2[:, :K] = t1[:, :K]
        t2[:, K: int(K + L / 2)] = (t1[:, K: : 2] + t1[:, K + 1: : 2] + 2 * np.abs(t1[:, K: : 2] - t1[:, K + 1: : 2])) / 3
        t3 = np.ones((N, M))
        for i in range(1, M):
            t3[:, i - 1] = r_sum(t2[:, list(range((i-1)*int(K/(M-1)), i*int(K/(M-1))))], np.ones((1, int(K/(M-1)))))
        t3[:, M - 1] = r_sum(t2[:, K: int(K + L / 2)], np.ones((1, int(L / 2))))
        x = np.zeros((N, M))
        for i in range(1, M):
            x[:, [i - 1]] = np.max([t3[:, [M - 1]], A[:, [i - 1]]], 0) * (t3[:, [i - 1]] - 0.5) + 0.5
        x[:, [M - 1]] = t3[:, [M - 1]]
        h = linear(x)
        f = np.tile(D * x[: ,[M - 1]], (1, M)) + np.tile(S, (N, 1)) * h
        pop.ObjV = f # 把求得的目标函数值赋值给种群pop的ObjV
    
    def calBest(self): # 计算全局最优解
        N = 10000 # 设置所要生成的全局最优解的个数
        X = np.hstack([np.array([np.linspace(0, 1, N)]).T, np.zeros((N, self.M - 2)) + 0.5, np.zeros((N, 1))])
        Point = linear(X)
        globalBestObjV = np.tile(np.array([list(range(2, 2 * self.M + 1, 2))]), (Point.shape[0], 1)) * Point
        return globalBestObjV

def linear(x):
    return np.fliplr(np.cumprod(np.hstack([np.ones((x.shape[0], 1)), x[:,:-1]]), 1)) * np.hstack([np.ones((x.shape[0], 1)), 1 - x[:, list(range(x.shape[1] - 1 - 1, -1, -1))]])

def s_linear(x, A):
    return np.abs(x - A) / np.abs(np.floor(A - x) + A)
    
def r_sum(x, w):
    Output = np.sum(x * np.tile(w, (x.shape[0], 1)), 1) / np.sum(w)
    return Output

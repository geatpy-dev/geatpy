# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ga

class tsp(ga.Problem): # 继承Problem父类
    def __init__(self, testName): # testName为测试集名称
        self.data=np.loadtxt("data/" + testName + ".csv",delimiter=",",usecols=(0,1)) # 读取城市坐标数据
        self.name = testName # 初始化name
        self.M = 1 # 初始化M（目标维数）
        self.maxormins = [1] * self.M # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = self.data.shape[0] # 初始化Dim（决策变量维数）
        self.varTypes = np.array([0] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [0] * self.Dim # 决策变量下界
        ub = [self.Dim - 1] * self.Dim # 决策变量上界
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
    
    def aimFuc(self, x, CV):
        x = x.copy()
        # 添加最后回到出发地
        X = np.hstack([x, x[:, [0]]])
        X = X.astype(int) 
        
        ObjV = [] # 存储所有种群个体对应的总路程
        for i in range(X.shape[0]):
            journey = self.data[X[i], :] # 按既定顺序到达的地点坐标
            distance = np.sum(np.sqrt(np.sum(np.diff(journey.T)**2, 0))) # 计算总路程
            ObjV.append(distance)
        ObjV = np.array([ObjV]).T
        
        return ObjV, CV

    
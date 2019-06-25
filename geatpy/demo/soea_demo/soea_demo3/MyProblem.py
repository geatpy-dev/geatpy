# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

"""
    一个带约束的单目标旅行商问题：
    有十座城市：A, B, C, D, E, F, G, H, I, J，坐标如下：
        X      Y
    [[0.4,  0.4439],
     [0.2439,0.1463],
     [0.1707,0.2293],
     [0.2293,0.761],
     [0.5171,0.9414],
     [0.8732,0.6536],
     [0.6878,0.5219],
     [0.8488,0.3609],
     [0.6683,0.2536],
     [0.6195,0.2634]]
    某旅行者从A城市出发，想逛遍所有城市，并且每座城市去且只去一次，最后要返回出发地，
而且需要从G地拿重要文件到D地，另外要从F地把公司的车开到E地，那么他应该如何设计行程方案，才能用
最短的路程来满足他的旅行需求？
    分析：在这个案例中，旅行者从A地出发，把其他城市走遍一次后回到A地，因此我们只需要考虑中间途
径的9个城市的访问顺序即可。这9个城市需要排列组合选出满足约束条件的最优的排列顺序作为最终的路线方案。
"""

class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self):
        self.name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        self.M = 1 # 初始化M（目标维数）
        self.maxormins = [1] # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = 9 # 初始化Dim（决策变量维数）
        self.varTypes = np.array([1] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [1] * self.Dim # 决策变量下界
        ub = [9] * self.Dim # 决策变量上界
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
        # 新增一个属性存储旅行地坐标
        self.place = np.array([[0.4,0.4439],
                               [0.2439,0.1463],
                               [0.1707,0.2293],
                               [0.2293,0.761],
                               [0.5171,0.9414],
                               [0.8732,0.6536],
                               [0.6878,0.5219],
                               [0.8488,0.3609],
                               [0.6683,0.2536],
                               [0.6195,0.2634]])
    
    def aimFuc(self, x, CV):
        # 添加从0地出发且最后回到出发地
        X = np.hstack([np.zeros((x.shape[0], 1)), x, np.zeros((x.shape[0], 1))]).astype(int)
        
        ObjV = [] # 存储所有种群个体对应的总路程
        for i in range(X.shape[0]):
            journey = self.place[X[i], :] # 按既定顺序到达的地点坐标
            distance = np.sum(np.sqrt(np.sum(np.diff(journey.T)**2, 0))) # 计算总路程
            ObjV.append(distance)
        ObjV = np.array([ObjV]).T
        # 找到违反约束条件的个体在种群中的索引，保存在向量exIdx中（如：若0、2、4号个体违反约束条件，则编程找出他们来）
        exIdx1 = np.where(np.where(x == 3)[1] - np.where(x == 6)[1] < 0)[0]
        exIdx2 = np.where(np.where(x == 4)[1] - np.where(x == 5)[1] < 0)[0]
        exIdx = np.unique(np.hstack([exIdx1, exIdx2]))
        CV[exIdx] = 1
        return [ObjV, CV]
    
    def calBest(self):
        realBestObjV = None
        return realBestObjV
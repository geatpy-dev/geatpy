# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea


class TestProblem(ea.Problem):  # 继承Problem父类
    def __init__(self, testName):  # testName为测试集名称
        name = testName  # 初始化name
        # 读取城市坐标数据
        self.places = np.loadtxt("data/" + testName + ".csv", delimiter=",", usecols=(0, 1))
        M = 1  # 初始化M（目标维数）
        Dim = self.places.shape[0]  # 初始化Dim（决策变量维数）
        maxormins = [1] * M  # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        varTypes = [0] * Dim  # 初始化varTypes（决策变量的类型，0：实数；1：整数）
        lb = [0] * Dim  # 决策变量下界
        ub = [Dim - 1] * Dim  # 决策变量上界
        lbin = [1] * Dim  # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1] * Dim  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)

    def aimFunc(self, pop):  # 目标函数
        x = pop.Phen.copy()  # 得到决策变量矩阵
        # 添加最后回到出发地
        X = np.hstack([x, x[:, [0]]]).astype(int)
        ObjV = []  # 存储所有种群个体对应的总路程
        for i in range(pop.sizes):
            journey = self.places[X[i], :]  # 按既定顺序到达的地点坐标
            distance = np.sum(np.sqrt(np.sum(np.diff(journey.T) ** 2, 0)))  # 计算总路程
            ObjV.append(distance)
        pop.ObjV = np.array([ObjV]).T

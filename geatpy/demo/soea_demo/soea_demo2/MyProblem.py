# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea
"""
该案例展示了一个带等式约束的连续型决策变量最大化目标的单目标优化问题。
该函数存在多个欺骗性很强的局部最优点。
max f = 4*x1 + 2*x2 + x3
s.t.
2*x1 + x2 - 1 <= 0
x1 + 2*x3 - 2 <= 0
x1 + x2 + x3 - 1 == 0
0 <= x1,x2 <= 1
0 < x3 < 2
"""
class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        M = 1 # 初始化M（目标维数）
        maxormins = [-1] # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 3 # 初始化Dim（决策变量维数）
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [0,0,0] # 决策变量下界
        ub = [1,1,2] # 决策变量上界
        lbin = [1,1,0] # 决策变量下边界
        ubin = [1,1,0] # 决策变量上边界
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        x1 = Vars[:, [0]]
        x2 = Vars[:, [1]]
        x3 = Vars[:, [2]]
        pop.ObjV = 4*x1 + 2*x2 + x3 # 计算目标函数值，赋值给pop种群对象的ObjV属性
        # 采用可行性法则处理约束
        pop.CV = np.hstack([2*x1 + x2 - 1,
                        x1 + 2*x3 - 2,
                        np.abs(x1 + x2 + x3 - 1)])
    
    def setBest(self): # 设置全局最优解
        globalBestObjV = np.array([[2.5]])
        return globalBestObjV
    
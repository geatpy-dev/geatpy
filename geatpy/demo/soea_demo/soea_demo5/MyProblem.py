# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

"""
该案例展示了一个需要混合编码种群来进化的最大化目标的单目标优化问题。
模型：
max f = sin(2x1) - cos(x2) + 2x3^2 -3x4 + (x5-3)^2 + 7x6
s.t.
-1.5 <= x1,x2 <= 2.5，
1 <= x3,x4,x5,x6 <= 7，且x3,x4,x5,x6为互不相等的整数。
分析：
该问题可以单纯用实整数编码'RI'来实现，但由于有一个”x3,x4,x5,x6互不相等“的约束，
因此把x3,x4,x5,x6用排列编码'P'，x1和x2采用实整数编码'RI'来求解会更好。
MyProblem是问题类，本质上是不需要管具体使用什么编码的，因此混合编码的设置在执行脚本main.py中进行而不是在此处。
"""

class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        M = 1 # 初始化M（目标维数）
        maxormins = [-1] # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 6 # 初始化Dim（决策变量维数）
        varTypes = [0,0,1,1,1,1] # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [-1.5,-1.5,1,1,1,1] # 决策变量下界
        ub = [2.5,2.5,7,7,7,7] # 决策变量上界
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
        pop.ObjV = np.sin(2*x1) - np.cos(x2) + 2*x3**2 -3*x4 + (x5-3)**2 + 7*x6 # 计算目标函数值，赋值给pop种群对象的ObjV属性
    
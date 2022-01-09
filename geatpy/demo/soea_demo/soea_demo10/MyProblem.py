# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

"""
    这是一个带等式约束的连续型决策变量最大化目标的单目标优化问题，
    该函数存在多个欺骗性很强的局部最优点。
    max f = 4*x1 + 2*x2 + x3
    s.t.
    2*x1 + x2 - 1 <= 0
    x1 + 2*x3 - 2 <= 0
    x1 + x2 + x3 - 1 == 0
    0 <= x1,x2 <= 1
    0 < x3 < 2
"""


class MyProblem(ea.Problem):  # 继承Problem父类
    def __init__(self):
        name = 'MyProblem'  # 初始化name（函数名称，可以随意设置）
        M = 1  # 初始化M（目标维数）
        maxormins = [-1]  # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 3  # 初始化Dim（决策变量维数）
        varTypes = [0] * Dim  # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [0, 0, 0]  # 决策变量下界
        ub = [1, 1, 2]  # 决策变量上界
        lbin = [1, 1, 0]  # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1, 1, 0]  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)

    def evalVars(self, Vars):  # 目标函数
        x1 = Vars[:, [0]]
        x2 = Vars[:, [1]]
        x3 = Vars[:, [2]]
        f = 4 * x1 + 2 * x2 + x3
        # 采用可行性法则处理约束
        CV = np.hstack([2 * x1 + x2 - 1,
                            x1 + 2 * x3 - 2,
                            np.abs(x1 + x2 + x3 - 1)])
        return f, CV

    def calReferObjV(self):  # 设定目标数参考值（本问题目标函数参考值设定为理论最优值）
        referenceObjV = np.array([[2.5]])
        return referenceObjV

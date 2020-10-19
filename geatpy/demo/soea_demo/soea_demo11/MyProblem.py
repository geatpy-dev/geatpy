# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

"""
该案例展示2个决策变量的单目标优化，决策变量的值将取自于一个设定好的变量集合。
max f = (-1+x1+((6-x2)*x2-2)*x2)**2+(-1+x1+((x2+2)*x2-10)*x2)**2
s.t.
x∈{1.1, 1, 0, 3, 5.5, 7.2, 9}
"""


class MyProblem(ea.Problem):  # 继承Problem父类
    def __init__(self):
        name = 'MyProblem'  # 初始化name（函数名称，可以随意设置）
        M = 1  # 初始化M（目标维数）
        maxormins = [1]  # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        self.var_set = np.array([1.1, 1, 0, 3, 5.5, 7.2, 9])  # 设定一个集合，要求决策变量的值取自于该集合
        Dim = 2  # 初始化Dim（决策变量维数）
        varTypes = [1] * Dim  # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [0, 0]  # 决策变量下界
        ub = [6, 6]  # 决策变量上界
        lbin = [1] * Dim  # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1] * Dim  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
        

    def aimFunc(self, pop):  # 目标函数
        Vars = pop.Phen.astype(np.int32)  # 强制类型转换确保元素是整数
        x1 = self.var_set[Vars[:, [0]]]  # 得到所有的x1组成的列向量
        x2 = self.var_set[Vars[:, [1]]]  # 得到所有的x2组成的列向量
        pop.ObjV = (-1 + x1 + ((6 - x2) * x2 - 2) * x2)**2 + (-1 + x1 + ((x2 + 2)*x2 - 10)*x2)**2  # 计算目标函数值，赋值给pop种群对象的ObjV属性

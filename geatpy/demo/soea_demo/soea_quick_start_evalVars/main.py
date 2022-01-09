# -*- coding: utf-8 -*-

import geatpy as ea
import numpy as np

"""

本案例展示了如何快速创建问题对象、快速调用算法求解一个含两个约束的单目标优化问题。
跟soea_quick_start_aimFunc案例不同的是，这里的目标函数定义成evalVars，而不是aimFunc。
区别：
    aimFunc(pop)传入的是一个种群对象。
    而evalVars(Vars)传入的是一个Numpy ndarray二维数组。其每一行的所有元素表示一组决策变量。
    函数evalVars可以有一个或两个返回值。第一个返回值表示目标函数值矩阵；第二个返回值表示违反约束程度矩阵。
    如果没有约束，则只需返回一个返回值即可。因此用法为：
    ObjV = evalVars(Vars)
    或
    ObjV, CV = evalVars(Vars)

"""

if __name__ == '__main__':
    # 构建问题
    r = 1  # 模拟该案例问题计算目标函数时需要用到的额外数据
    def evalVars(Vars):  # 定义目标函数（含约束）
        ObjV = np.sum((Vars - r) ** 2, 1, keepdims=True)  # 计算目标函数值
        x1 = Vars[:, [0]]  # 把Vars的第0列取出来
        x2 = Vars[:, [1]]  # 把Vars的第1列取出来
        CV = np.hstack([(x1 - 0.5) ** 2 - 0.25,
                        (x2 - 1) ** 2 - 1])  # 计算违反约束程度
        return ObjV, CV  # 返回目标函数值矩阵和违反约束程度矩阵

    problem = ea.Problem(name='soea quick start demo',
                         M=1,  # 目标维数
                         maxormins=[1],  # 目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标
                         Dim=5,  # 决策变量维数
                         varTypes=[0, 0, 1, 1, 1],  # 决策变量的类型列表，0：实数；1：整数
                         lb=[-1, 1, 2, 1, 0],  # 决策变量下界
                         ub=[1, 4, 5, 2, 1],  # 决策变量上界
                         evalVars=evalVars)
    # 构建算法
    algorithm = ea.soea_SEGA_templet(problem,
                                     ea.Population(Encoding='RI', NIND=20),
                                     MAXGEN=50,  # 最大进化代数。
                                     logTras=1,  # 表示每隔多少代记录一次日志信息，0表示不记录。
                                     trappedValue=1e-6,  # 单目标优化陷入停滞的判断阈值。
                                     maxTrappedCount=10)  # 进化停滞计数器最大上限值。
    # 求解
    res = ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=True, drawLog=False, saveFlag=True)
    print(res)

# -*- coding: utf-8 -*-

import geatpy as ea
import numpy as np

"""

本案例展示了如何快速创建问题对象、快速调用算法求解一个含两个约束的单目标优化问题。
与soea_quick_start_evalVars案例不同的是，本案例采用循环计算每个个体的目标函数值和违反约束程度值。

"""

if __name__ == '__main__':
    # 构建问题
    r = 1  # 目标函数需要用到的额外数据
    def evalVars(Vars):  # 定义目标函数（含约束）
        f = []  # 存储目标函数值
        CV = []  # 存储违反约束程度
        for i in range(Vars.shape[0]):  # 遍历每一组决策变量，计算对应的目标函数值
            f.append(np.sum((Vars[i, :] - r) ** 2))  # 用Vars[i, :] 取出每一组决策变量
            x1 = Vars[i, 0]
            x2 = Vars[i, 1]
            CV.append(np.array([(x1 - 0.5)**2 - 0.25,
                                (x2 - 1)**2 - 1]))
        return np.vstack(f), np.vstack(CV)  # 返回目标函数值矩阵和违反约束程度矩阵

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

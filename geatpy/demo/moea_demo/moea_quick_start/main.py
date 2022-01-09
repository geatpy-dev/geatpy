# -*- coding: utf-8 -*-

import geatpy as ea
import numpy as np

"""

本案例展示了如何快速创建多目标优化问题对象、快速调用算法求解。

"""

if __name__ == '__main__':
    # 构建问题
    def evalVars(Vars):  # 定义目标函数（含约束）
        f1 = Vars ** 2
        f2 = (Vars - 2) ** 2
        ObjV = np.hstack([f1, f2])
        CV = -Vars ** 2 + 2.5 * Vars - 1.5
        return ObjV, CV


    problem = ea.Problem(name='moea quick start',
                         M=2,
                         maxormins=[1, 1],
                         Dim=1,
                         varTypes=[0],
                         lb=[-10],
                         ub=[10],
                         evalVars=evalVars)
    # 构建算法
    algorithm = ea.moea_NSGA2_templet(problem,
                                      ea.Population(Encoding='RI', NIND=20),
                                      MAXGEN=50,  # 最大进化代数
                                      logTras=1)  # 表示每隔多少代记录一次日志信息，0表示不记录。
    # 求解
    res = ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=True, drawLog=True, saveFlag=True)
    print(res)

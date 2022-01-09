# -*- coding: utf-8 -*-
import numpy as np
from MyProblem import MyProblem  # 导入自定义问题接口
import geatpy as ea  # import geatpy

"""
    该案例展示了一个带等式约束的连续型决策变量最大化目标的单目标优化问题。问题的定义详见MyProblem.py。
    与soea_demo2不同之处是采用先验知识帮助进化优化。
"""

if __name__ == '__main__':
    # 实例化问题对象
    problem = MyProblem()  # 生成问题对象
    # 快速构建算法
    algorithm = ea.soea_DE_currentToBest_1_bin_templet(problem,
                                                       ea.Population(Encoding='RI', NIND=20),
                                                       MAXGEN=400,  # 最大进化代数。
                                                       logTras=0)  # 表示每隔多少代记录一次日志信息，0表示不记录。
    algorithm.mutOper.F = 0.7  # 差分进化中的参数F。
    algorithm.recOper.XOVR = 0.7  # 交叉概率。
    # 先验知识
    prophetVars = np.array([[0.4, 0.2, 0.4]])  # 假设已知[0.4, 0.2, 0.4]为一组比较优秀的变量。
    # 求解
    res = ea.optimize(algorithm, prophet=prophetVars, verbose=True, drawing=1, outputMsg=True, drawLog=True, saveFlag=True)
    print(res)

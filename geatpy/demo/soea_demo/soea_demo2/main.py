# -*- coding: utf-8 -*-
from MyProblem import MyProblem  # 导入自定义问题接口
import geatpy as ea  # import geatpy

"""
    该案例展示了一个带等式约束的连续型决策变量最大化目标的单目标优化问题的求解。问题的定义详见MyProblem.py。
"""

if __name__ == '__main__':
    # 实例化问题对象
    problem = MyProblem()
    # 构建算法
    algorithm = ea.soea_DE_rand_1_bin_templet(problem,
                                              ea.Population(Encoding='RI', NIND=100),
                                              MAXGEN=500,  # 最大进化代数。
                                              logTras=1)  # 表示每隔多少代记录一次日志信息，0表示不记录。
    algorithm.mutOper.F = 0.5  # 差分进化中的参数F
    algorithm.recOper.XOVR = 0.7  # 重组概率
    # 求解
    res = ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=True, drawLog=False, saveFlag=True)
    print(res)

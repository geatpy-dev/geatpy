# -*- coding: utf-8 -*-
from MyProblem import MyProblem  # 导入自定义问题接口
import geatpy as ea  # import geatpy

"""
该案例展示了一个简单的连续型决策变量最大化目标的单目标优化问题的求解。问题的定义详见MyProblem.py。
"""

if __name__ == '__main__':
    # 实例化问题对象
    problem = MyProblem()
    # 构建算法
    algorithm = ea.soea_SEGA_templet(problem,
                                     ea.Population(Encoding='RI', NIND=40),
                                     MAXGEN=25,  # 最大进化代数。
                                     logTras=1,  # 表示每隔多少代记录一次日志信息，0表示不记录。
                                     trappedValue=1e-6,  # 单目标优化陷入停滞的判断阈值。
                                     maxTrappedCount=10)  # 进化停滞计数器最大上限值。
    # 求解
    res = ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=True, drawLog=False, saveFlag=True)
    print(res)

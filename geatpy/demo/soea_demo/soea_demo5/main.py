# -*- coding: utf-8 -*-
from MyProblem import MyProblem  # 导入自定义问题接口
import geatpy as ea  # import geatpy

"""
该案例展示了一个需要混合编码种群来进化的最大化目标的单目标优化问题。问题的定义详见MyProblem.py。
分析：
该问题可以单纯用实整数编码'RI'来实现，但由于有一个”x3,x4,x5,x6互不相等“的约束，
因此把x3,x4,x5,x6用排列编码'P'，x1和x2采用实整数编码'RI'来求解会更好。
MyProblem是问题类，本质上是不需要管具体使用什么编码的，因此混合编码的设置在执行脚本main.py中进行而不是在此处。
"""

if __name__ == '__main__':
    # 实例化问题对象
    problem = MyProblem()
    # 快速构建算法
    algorithm = ea.soea_psy_EGA_templet(problem,
                                        ea.PsyPopulation(Encodings=['RI', 'P'], NIND=40, EncoIdxs=[[0, 1], [2, 3, 4, 5]]),
                                        MAXGEN=25,  # 最大进化代数
                                        logTras=1)  # 表示每隔多少代记录一次日志信息，0表示不记录。
    # 求解
    res = ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=True, drawLog=False, saveFlag=True)
    print(res)

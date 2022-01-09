# -*- coding: utf-8 -*-
from MyProblem import MyProblem  # 导入自定义问题接口
import geatpy as ea  # import geatpy

"""
该案例是soea_demo1的拓展，展示了如何利用多种群来进行单目标优化。问题的定义详见MyProblem.py。
在执行脚本main.py中，通过调用带"multi"字样的进化算法类来进行多种群进化优化。
"""

if __name__ == '__main__':
    # 实例化问题对象
    problem = MyProblem()
    # 种群设置
    Encoding = 'RI'  # 编码方式
    NINDs = [5, 10, 15, 20]  # 种群规模
    population = []  # 创建种群列表
    for i in range(len(NINDs)):
        Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders)  # 创建区域描述器。
        population.append(ea.Population(Encoding, Field, NINDs[i]))  # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）。
    # 构建算法
    algorithm = ea.soea_multi_SEGA_templet(problem,
                                           population,
                                           MAXGEN=30,  # 最大进化代数。
                                           logTras=1,  # 表示每隔多少代记录一次日志信息，0表示不记录。
                                           trappedValue=1e-6,  # 单目标优化陷入停滞的判断阈值。
                                           maxTrappedCount=5)  # 进化停滞计数器最大上限值。
    # 求解
    res = ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=True, drawLog=False, saveFlag=False)
    print(res)

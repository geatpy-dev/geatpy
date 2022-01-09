# -*- coding: utf-8 -*-
from MyProblem import MyProblem  # 导入自定义问题接口
import geatpy as ea  # import geatpy

"""
描述:
    该案例是moea_demo1的另一个版本，展示了如何定义aimFunc()而不是evalVars()来计算目标函数和违反约束程度值。【见MyProblem.py】
    同时展示如何定义outFunc()，用于让算法在每一次进化时调用该outFunc()函数。
"""

if __name__ == '__main__':
    #  实例化问题对象
    problem = MyProblem()
    # 定义outFunc()函数
    def outFunc(alg, pop):  # alg 和 pop为outFunc的固定输入参数，分别为算法对象和每次迭代的种群对象。
        print('第 %d 代' % alg.currentGen)
    # 构建算法
    algorithm = ea.moea_NSGA2_templet(problem,
                                      ea.Population(Encoding='RI', NIND=50),
                                      MAXGEN=200,  # 最大进化代数
                                      logTras=1,  # 表示每隔多少代记录一次日志信息，0表示不记录。
                                      outFunc=outFunc)
    # 求解
    res = ea.optimize(algorithm, verbose=False, drawing=1, outputMsg=True, drawLog=True, saveFlag=False)
    print(res)

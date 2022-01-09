# -*- coding: utf-8 -*-
from MyProblem import MyProblem  # 导入自定义问题接口
import geatpy as ea  # import geatpy

if __name__ == '__main__':
    # 实例化问题对象
    problem = MyProblem()
    # 构建算法
    algorithm = ea.moea_psy_NSGA2_templet(problem,
                                          ea.PsyPopulation(Encodings=['RI', 'P'], NIND=30, EncoIdxs=[[0], [1,2,3,4]]),
                                          MAXGEN=200,  # 最大进化代数
                                          logTras=0)  # 表示每隔多少代记录一次日志信息，0表示不记录。
    # 求解
    res = ea.optimize(algorithm, verbose=False, drawing=1, outputMsg=True, drawLog=False, saveFlag=True)
    print(res)

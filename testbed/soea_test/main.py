# -*- coding: utf-8 -*-
import geatpy as ea  # import geatpy

if __name__ == '__main__':
    # 问题对象
    problem = ea.benchmarks.Ackley(30)
    # 构建算法
    algorithm = ea.soea_DE_rand_1_bin_templet(
        problem,
        ea.Population(Encoding='RI', NIND=20),
        MAXGEN=1000,  # 最大进化代数。
        logTras=1)  # 表示每隔多少代记录一次日志信息，0表示不记录。
    algorithm.mutOper.F = 0.5  # 差分进化中的参数F
    algorithm.recOper.XOVR = 0.2  # 差分进化中的参数Cr
    # 求解
    res = ea.optimize(algorithm,
                      verbose=True,
                      drawing=1,
                      outputMsg=True,
                      drawLog=True,
                      saveFlag=True,
                      dirName='result')
    print(res)

# -*- coding: utf-8 -*-
from MyProblem import MyProblem  # 导入自定义问题接口
import geatpy as ea  # import geatpy

"""
该案例展示了如何利用进化算法+多进程/多线程来搜索SVM中的两个参数：C和Gamma的最佳值。问题的定义详见MyProblem.py。
在执行脚本main.py中设置PoolType字符串来控制采用的是多进程还是多线程。
注意：使用多进程时，程序必须以“if __name__ == '__main__':”作为入口，
      这个是multiprocessing的多进程模块的硬性要求。
"""

if __name__ == '__main__':
    # 实例化问题对象
    problem = MyProblem(PoolType='Thread')  # 设置采用多线程，若修改为: PoolType = 'Process'，则表示用多进程
    # 构建算法
    algorithm = ea.soea_DE_rand_1_bin_templet(problem,
                                              ea.Population(Encoding='RI', NIND=50),
                                              MAXGEN=30,  # 最大进化代数。
                                              logTras=1,  # 表示每隔多少代记录一次日志信息，0表示不记录。
                                              trappedValue=1e-6,  # 单目标优化陷入停滞的判断阈值。
                                              maxTrappedCount=10)  # 进化停滞计数器最大上限值。
    # 求解
    res = ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=True, drawLog=False, saveFlag=True)
    # 检验结果
    if res['success']:
        problem.test(C=res['Vars'][0, 0], G=res['Vars'][0, 1])

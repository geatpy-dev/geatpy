# -*- coding: utf-8 -*-
from MyProblem import MyProblem  # 导入自定义问题接口
import geatpy as ea  # import geatpy

"""
该案例展示了一个利用单目标进化算法实现句子匹配的应用实例。问题的定义详见MyProblem.py。
"""

if __name__ == '__main__':
    # 实例化问题对象
    problem = MyProblem()
    # 快速构建算法
    algorithm = ea.soea_DE_rand_1_L_templet(problem,
                                            ea.Population(Encoding='RI', NIND=50),
                                            MAXGEN=2000,  # 最大进化代数
                                            logTras=1)  # 表示每隔多少代记录一次日志信息，0表示不记录。
    # 求解
    res = ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=False, drawLog=False, saveFlag=True)
    print('最佳目标函数值：%s' % res['ObjV'][0][0])
    print('搜索到的句子为：')
    for num in res['Vars'][0, :]:
        print(chr(int(num)), end='')

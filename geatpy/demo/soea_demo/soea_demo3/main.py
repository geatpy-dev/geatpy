# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from MyProblem import MyProblem  # 导入自定义问题接口
import geatpy as ea  # import geatpy

"""
    该案例展示了一个带约束的单目标旅行商问题的求解。问题的定义详见MyProblem.py。
"""

if __name__ == '__main__':
    # 实例化问题对象
    problem = MyProblem()
    # 构建算法
    algorithm = ea.soea_SEGA_templet(problem,
                                     ea.Population(Encoding='P', NIND=50),
                                     MAXGEN=200,  # 最大进化代数
                                     logTras=1)  # 表示每隔多少代记录一次日志信息，0表示不记录。
    algorithm.mutOper.Pm = 0.5  # 变异概率
    # 求解
    res = ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=True, drawLog=False, saveFlag=True)
    # 绘制路线图
    if res['success']:
        print('最短路程为：%s' % res['ObjV'][0][0])
        print('最佳路线为：')
        best_journey = np.hstack([0, res['Vars'][0, :], 0])
        for i in range(len(best_journey)):
            print(int(best_journey[i]), end=' ')
        print()
        # 绘图
        plt.figure()
        plt.plot(problem.places[best_journey.astype(int), 0], problem.places[best_journey.astype(int), 1], c='black')
        plt.plot(problem.places[best_journey.astype(int), 0], problem.places[best_journey.astype(int), 1], 'o',
                 c='black')
        for i in range(len(best_journey)):
            plt.text(problem.places[int(best_journey[i]), 0], problem.places[int(best_journey[i]), 1],
                     chr(int(best_journey[i]) + 65), fontsize=20)
        plt.grid(True)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig('roadmap.svg', dpi=600, bbox_inches='tight')
        plt.show()
    else:
        print('没找到可行解。')

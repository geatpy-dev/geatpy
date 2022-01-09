# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import geatpy as ea  # import geatpy

if __name__ == '__main__':
    # 实例化问题对象
    problem = ea.benchmarks.TSP('att48')
    # 构建算法
    algorithm = ea.soea_studGA_templet(problem,
                                      ea.Population(Encoding='P', NIND=100),
                                      MAXGEN=1000,  # 最大进化代数
                                      logTras=1)  # 表示每隔多少代记录一次日志信息
    algorithm.mutOper.Pm = 0.5  # 变异概率
    # 求解
    saveDirName = 'result'
    res = ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=True, drawLog=False, saveFlag=True, dirName=saveDirName)
    # 绘制路线图
    if res['success']:
        print('最短路程为：%s' % res['ObjV'][0][0])
        print('最佳路线为：')
        best_journey = np.hstack([res['Vars'][0, :], res['Vars'][0, 0]])
        for i in range(len(best_journey)):
            print(int(best_journey[i]), end=' ')
        print()
        # 绘图
        plt.figure()
        plt.plot(problem.places[best_journey.astype(int), 0], problem.places[best_journey.astype(int), 1], c='black')
        plt.plot(problem.places[best_journey.astype(int), 0], problem.places[best_journey.astype(int), 1], 'o', c='black')
        for i in range(len(best_journey)):
            plt.text(problem.places[int(best_journey[i]), 0], problem.places[int(best_journey[i]), 1],
                     int(best_journey[i]), fontsize=15)
        plt.grid(True)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig(saveDirName + '/roadmap.svg', dpi=600, bbox_inches='tight')
        plt.show()
    else:
        print('没找到可行解。')

# -*- coding: utf-8 -*-
import geatpy as ea  # import geatpy
import numpy as np
import matplotlib.pyplot as plt
from tsp import TestProblem

if __name__ == '__main__':
    """===============================实例化问题对象============================"""
    problem = TestProblem('att48')  # 生成问题对象
    """=================================种群设置==============================="""
    Encoding = 'P'  # 编码方式
    NIND = 100  # 种群规模
    Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders)  # 创建区域描述器
    population = ea.Population(Encoding, Field, NIND)  # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
    """===============================算法参数设置=============================="""
    myAlgorithm = ea.soea_studGA_templet(problem, population)  # 实例化一个算法模板对象
    myAlgorithm.MAXGEN = 1000  # 最大进化代数
    myAlgorithm.logTras = 1  # 设置每隔多少代记录日志，若设置成0则表示不记录日志
    myAlgorithm.verbose = True  # 设置是否打印输出日志信息
    myAlgorithm.drawing = 1  # 设置绘图方式（0：不绘图；1：绘制结果图；2：绘制目标空间过程动画；3：绘制决策空间过程动画）
    """==========================调用算法模板进行种群进化========================="""
    [BestIndi, population] = myAlgorithm.run()  # 执行算法模板，得到最优个体以及最后一代种群
    BestIndi.save()  # 把最优个体的信息保存到文件中
    """=================================输出结果==============================="""
    print('评价次数：%s' % myAlgorithm.evalsNum)
    print('时间已过 %s 秒' % myAlgorithm.passTime)
    if BestIndi.sizes != 0:
        print('最短路程为：%s' % BestIndi.ObjV[0][0])
        print('最佳路线为：')
        best_journey = np.hstack([BestIndi.Phen[0, :], BestIndi.Phen[0, 0]])
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
                     int(best_journey[i]), fontsize=15)
        plt.grid(True)
        plt.xlabel('x坐标')
        plt.ylabel('y坐标')
        plt.savefig('roadmap.svg', dpi=600, bbox_inches='tight')
    else:
        print('没找到可行解。')

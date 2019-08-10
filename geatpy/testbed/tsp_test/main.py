# -*- coding: utf-8 -*-
import geatpy as ea # import geatpy
import numpy as np
import matplotlib.pyplot as plt
from tsp import TestProblem

"""===============================实例化问题对象=============================="""
problem = TestProblem('att48')        # 生成问题对象
"""=================================种群设置================================="""
Encoding = 'P'                        # 编码方式
NIND = 100                            # 种群规模
Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders) # 创建区域描述器
population = ea.Population(Encoding, Field, NIND) # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
"""===============================算法参数设置==============================="""
myAlgorithm = ea.soea_SEGA_templet(problem, population) # 实例化一个算法模板对象
myAlgorithm.MAXGEN = 2000             # 最大进化代数
myAlgorithm.drawing = 1               # 设置绘图方式（0：不绘图；1：绘制结果图；2：绘制过程动画）
"""==========================调用算法模板进行种群进化=========================="""
[population, obj_trace, var_trace] = myAlgorithm.run() # 执行算法模板，得到最后一代种群以及进化记录器
population.save()                     # 把最后一代种群的信息保存到文件中
# 输出结果
best_gen = np.argmin(obj_trace[:, 1]) # 记录最优种群是在哪一代
best_ObjV = np.min(obj_trace[:, 1])
print('最短路程为：%s'%(best_ObjV))
print('最佳路线为：')
best_journey = np.hstack([var_trace[best_gen, :], var_trace[best_gen, 0]])
for i in range(len(best_journey)):
    print(int(best_journey[i]), end = ' ')
print()
print('有效进化代数：%s'%(obj_trace.shape[0]))
print('最优的一代是第 %s 代'%(best_gen + 1))
print('评价次数：%s'%(myAlgorithm.evalsNum))
print('时间已过 %s 秒'%(myAlgorithm.passTime))
# 绘图
plt.figure()
plt.plot(problem.data[best_journey.astype(int), 0], problem.data[best_journey.astype(int), 1], c = 'black')
for i in range(len(best_journey)):
    plt.text(problem.data[int(best_journey[i]), 0], problem.data[int(best_journey[i]), 1], int(best_journey[i]), fontsize=15)
plt.xlabel('x坐标')
plt.ylabel('y坐标')
plt.savefig('roadmap.svg', dpi=600, bbox_inches='tight')

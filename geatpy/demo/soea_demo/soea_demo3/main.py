# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea # import geatpy
import matplotlib.pyplot as plt
from MyProblem import MyProblem # 导入自定义问题接口

"""==================================实例化问题对象================================"""
problem = MyProblem() # 生成问题对象
"""==================================种群设置================================"""
Encoding = 'P'             # 编码方式，采用排列编码
conordis = 1               # 表示染色体解码后得到的变量是离散的
NIND = 50                  # 种群规模
Field = ea.crtfld(Encoding, conordis, problem.ranges, problem.borders) # 创建区域描述器
population = ea.Population(Encoding, conordis, Field, NIND) # 实例化种群对象（此时种群还没被真正初始化）
"""==================================算法参数设置================================"""
myAlgorithm = ea.soea_SEGA_templet(problem, population) # 实例化一个算法模板对象
myAlgorithm.MAXGEN = 100 # 最大遗传代数
myAlgorithm.drawing = 1
"""=======================调用算法模板进行种群进化=============================="""
[population, obj_trace, var_trace] = myAlgorithm.run() # 执行算法模板
"""===============================输出结果及绘图=============================="""
# 输出结果
best_gen = np.argmin(obj_trace[:, 1]) # 记录最优种群是在哪一代
best_ObjV = np.min(obj_trace[:, 1])
print('最短路程为：%s'%(best_ObjV))
print('最佳路线为：')
best_journey = np.hstack([0, var_trace[best_gen, :], 0])
for i in range(len(best_journey)):
    print(int(best_journey[i]), end = ' ')
print()
print('有效进化代数：%s'%(obj_trace.shape[0]))
print('最优的一代是第 %s 代'%(best_gen + 1))
print('评价次数：%s'%(myAlgorithm.evalsNum))
print('时间已过 %s 秒'%(myAlgorithm.passTime))
# 绘图
plt.plot(problem.place[best_journey.astype(int), 0], problem.place[best_journey.astype(int), 1], c = 'black')
for i in range(len(best_journey) - 1):
    plt.annotate('', xy=(problem.place[int(best_journey[i+1]), 0], problem.place[int(best_journey[i+1]), 1]), xytext=(problem.place[int(best_journey[i]), 0], problem.place[int(best_journey[i]), 1]),arrowprops=dict(facecolor='black', width=1, headwidth = 5, shrink=0.01))
for i in range(len(best_journey)):
    plt.text(problem.place[int(best_journey[i]), 0], problem.place[int(best_journey[i]), 1], chr(int(best_journey[i]) + 65), fontsize=20)
plt.xlabel('x坐标')
plt.ylabel('y坐标')

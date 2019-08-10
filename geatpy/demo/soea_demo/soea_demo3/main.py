# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea # import geatpy
import matplotlib.pyplot as plt
from MyProblem import MyProblem # 导入自定义问题接口

"""==================================实例化问题对象================================"""
problem = MyProblem() # 生成问题对象
"""==================================种群设置================================"""
Encoding = 'P'        # 编码方式，采用排列编码
NIND = 50             # 种群规模
Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders) # 创建区域描述器
population = ea.Population(Encoding, Field, NIND) # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
"""==================================算法参数设置================================"""
myAlgorithm = ea.soea_SEGA_templet(problem, population) # 实例化一个算法模板对象
myAlgorithm.MAXGEN = 200 # 最大进化代数
myAlgorithm.drawing = 1
myAlgorithm.pm = 0.5
"""=======================调用算法模板进行种群进化=============================="""
[population, obj_trace, var_trace] = myAlgorithm.run() # 执行算法模板
population.save() # 把最后一代种群的信息保存到文件中
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
plt.figure()
plt.plot(problem.place[best_journey.astype(int), 0], problem.place[best_journey.astype(int), 1], c = 'black')
for i in range(len(best_journey) - 1):
    plt.annotate('', xy=(problem.place[int(best_journey[i+1]), 0], problem.place[int(best_journey[i+1]), 1]), xytext=(problem.place[int(best_journey[i]), 0], problem.place[int(best_journey[i]), 1]),arrowprops=dict(facecolor='black', width=1, headwidth = 5, shrink=0.01))
for i in range(len(best_journey)):
    plt.text(problem.place[int(best_journey[i]), 0], problem.place[int(best_journey[i]), 1], chr(int(best_journey[i]) + 65), fontsize=20)
plt.xlabel('x坐标')
plt.ylabel('y坐标')
plt.savefig('roadmap.svg', dpi=600, bbox_inches='tight')
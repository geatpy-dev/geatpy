# -*- coding: utf-8 -*-
import geatpy as ea # import geatpy
import numpy as np

"""==================================实例化问题对象================================"""
problemName = 'Shubert' # 问题名称
fileName = problemName  # 这里因为目标函数写在与之同名的文件里，所以文件名也是问题名称
MyProblem = getattr(__import__(fileName), problemName) # 获得自定义问题类
problem = MyProblem()   # 生成问题对象
"""==================================种群设置================================"""
Encoding = 'RI'         # 编码方式
NIND = 20               # 种群规模
precisions = [50] * problem.Dim # 编码精度（适用于二进制/格雷编码）
Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders, precisions) # 创建区域描述器
population = ea.Population(Encoding, Field, NIND) # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
"""==================================算法参数设置================================"""
myAlgorithm = ea.soea_DE_rand_1_bin_templet(problem, population) # 实例化一个算法模板对象
myAlgorithm.MAXGEN = 1000 # 最大遗传代数
myAlgorithm.F = 0.5
myAlgorithm.pc = 0.2
myAlgorithm.drawing = 1 # 设置绘图方式（0：不绘图；1：绘制结果图；2：绘制过程动画）
"""=======================调用算法模板进行种群进化=============================="""
[population, obj_trace, var_trace] = myAlgorithm.run() # 执行算法模板，得到最后一代种群以及进化记录器
population.save() # 把最后一代种群的信息保存到文件中
# 输出结果
best_gen = np.argmin(obj_trace[:, 1]) # 记录最优种群是在哪一代
best_ObjV = np.min(obj_trace[:, 1])
print('最优的目标函数值为：%s'%(best_ObjV))
print('最优的控制变量值为：')
for i in range(var_trace.shape[1]):
    print(var_trace[best_gen, i])
print('有效进化代数：%s'%(obj_trace.shape[0]))
print('最优的一代是第 %s 代'%(best_gen + 1))
print('评价次数：%s'%(myAlgorithm.evalsNum))
print('时间已过 %s 秒'%(myAlgorithm.passTime))

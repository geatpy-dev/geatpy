# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea # import geatpy
from MyProblem import MyProblem # 导入自定义问题接口

"""================================实例化问题对象============================="""
problem = MyProblem() # 生成问题对象
"""==================================种群设置================================"""
NIND = 40             # 种群规模
# 创建区域描述器，这里需要创建两个，前2个变量用RI编码，剩余变量用排列编码
Encodings = ['RI', 'P']
Field1 = ea.crtfld(Encodings[0], problem.varTypes[:2], problem.ranges[:,:2], problem.borders[:,:2])
Field2 = ea.crtfld(Encodings[1], problem.varTypes[2:], problem.ranges[:,2:], problem.borders[:,2:])
Fields = [Field1, Field2]
population = ea.PsyPopulation(Encodings, Fields, NIND) # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
"""================================算法参数设置=============================="""
myAlgorithm = ea.soea_psy_SEGA_templet(problem, population) # 实例化一个算法模板对象
myAlgorithm.MAXGEN = 25 # 最大进化代数
"""===========================调用算法模板进行种群进化========================="""
[population, obj_trace, var_trace] = myAlgorithm.run() # 执行算法模板
population.save() # 把最后一代种群的信息保存到文件中
# 输出结果
best_gen = np.argmax(obj_trace[:, 1]) # 记录最优种群是在哪一代
best_ObjV = obj_trace[best_gen, 1]
print('最优的目标函数值为：%s'%(best_ObjV))
print('最优的控制变量值为：')
for i in range(var_trace.shape[1]):
    print(var_trace[best_gen, i])
print('有效进化代数：%s'%(obj_trace.shape[0]))
print('最优的一代是第 %s 代'%(best_gen + 1))
print('评价次数：%s'%(myAlgorithm.evalsNum))
print('时间已过 %s 秒'%(myAlgorithm.passTime))

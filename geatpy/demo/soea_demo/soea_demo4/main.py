# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea # import geatpy
from MyProblem import MyProblem # 导入自定义问题接口

"""==================================实例化问题对象================================"""
problem = MyProblem() # 生成问题对象
"""==================================种群设置================================"""
Encoding = 'RI'       # 编码方式
NIND = 50             # 种群规模
Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders) # 创建区域描述器
population = ea.Population(Encoding, Field, NIND) # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
"""==================================算法参数设置================================"""
myAlgorithm = ea.soea_DE_rand_1_L_templet(problem, population) # 实例化一个算法模板对象
myAlgorithm.MAXGEN = 2000  # 最大进化代数
"""=======================调用算法模板进行种群进化=============================="""
[population, obj_trace, var_trace] = myAlgorithm.run() # 执行算法模板
population.save() # 把最后一代种群的信息保存到文件中
# 输出结果
best_gen = np.argmin(obj_trace[:, 1]) # 记录最优种群是在哪一代
best_ObjV = obj_trace[best_gen, 1]
print('最优的目标函数值为：%s'%(best_ObjV))
print('有效进化代数：%s'%(obj_trace.shape[0]))
print('最优的一代是第 %s 代'%(best_gen + 1))
print('评价次数：%s'%(myAlgorithm.evalsNum))
print('时间已过 %s 秒'%(myAlgorithm.passTime))
for num in var_trace[best_gen, :]:
    print(chr(int(num)), end = '')

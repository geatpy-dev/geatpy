# -*- coding: utf-8 -*-
import geatpy as ea # import geatpy
from MyProblem import MyProblem # 导入自定义问题接口

"""==================================实例化问题对象================================"""
problem = MyProblem()    # 生成问题对象
"""==================================种群设置================================"""
NIND = 30                # 种群规模
# 创建区域描述器，这里需要创建两个，前1个变量用RI编码，剩余变量用排列编码
Encodings = ['RI', 'P']  # 编码方式列表
Field1 = ea.crtfld(Encodings[0], problem.varTypes[:1], problem.ranges[:,:1], problem.borders[:,:1]) # 创建区域描述器
Field2 = ea.crtfld(Encodings[1], problem.varTypes[1:], problem.ranges[:,1:], problem.borders[:,1:])
Fields = [Field1, Field2]
population = ea.PsyPopulation(Encodings, Fields, NIND) # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
"""==================================算法参数设置================================"""
myAlgorithm = ea.moea_psy_NSGA2_templet(problem, population) # 实例化一个算法模板对象
myAlgorithm.MAXGEN = 200    # 最大进化代数
"""=======================调用算法模板进行种群进化=============================="""
NDSet = myAlgorithm.run()   # 执行算法模板，得到帕累托最优解集NDSet
NDSet.save()                # 把结果保存到文件中
# 输出
print('用时：%s 秒'%(myAlgorithm.passTime))
print('非支配个体数：%s 个'%(NDSet.sizes))
print('单位时间找到帕累托前沿点个数：%s 个'%(int(NDSet.sizes // myAlgorithm.passTime)))

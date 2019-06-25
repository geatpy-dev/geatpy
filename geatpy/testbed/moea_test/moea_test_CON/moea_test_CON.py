# -*- coding: utf-8 -*-
import geatpy as ea # import geatpy

"""==================================实例化问题对象================================"""
problemName = 'CON' # 目标函数名
fileName = 'MyProblem' # 目标函数所在的文件名
MyProblem = getattr(__import__(fileName), problemName) # 获得自定义问题类
problem = MyProblem() # 生成问题对象

"""==================================种群设置================================"""
Encoding = 'R'             # 编码方式
conordis = 0               # 表示染色体解码后得到的变量是连续的
NIND = 200                 # 种群规模
FieldDR = ea.crtfld(Encoding, conordis, problem.ranges, problem.borders) # 创建区域描述器
population = ea.Population(Encoding, conordis, FieldDR, NIND) # 实例化种群对象（此时种群还没被真正初始化）

"""==================================算法参数设置================================"""
myAlgorithm = ea.moea_NSGA2_templet(problem, population) # 实例化一个算法模板对象
myAlgorithm.MAXGEN = 50 # 最大遗传代数
"""=======================调用算法模板进行种群进化=============================="""
NDSet = myAlgorithm.run() # 执行算法模板，得到帕累托最优解集NDSet
# 输出
print('用时：%s 秒'%(myAlgorithm.passTime))
print('非支配个体数：%s 个'%(NDSet.sizes))
print('单位时间找到帕累托前沿点个数：%s 个'%(int(NDSet.sizes // myAlgorithm.passTime)))

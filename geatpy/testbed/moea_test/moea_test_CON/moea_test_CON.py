# -*- coding: utf-8 -*-
import geatpy as ea # import geatpy

"""==================================实例化问题对象================================"""
problemName = 'CON'       # 问题名称
fileName = 'MyProblem'    # 问题类所在的文件名
MyProblem = getattr(__import__(fileName), problemName) # 获得自定义问题类
problem = MyProblem()     # 生成问题对象
"""==================================种群设置================================"""
Encoding = 'RI'           # 编码方式
NIND = 200                # 种群规模
Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders) # 创建区域描述器
population = ea.Population(Encoding, Field, NIND) # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
"""==================================算法参数设置================================"""
myAlgorithm = ea.moea_NSGA2_templet(problem, population) # 实例化一个算法模板对象
myAlgorithm.MAXGEN = 50   # 最大进化代数
myAlgorithm.drawing = 1   # 设置绘图方式（0：不绘图；1：绘制结果图；2：绘制过程动画）
"""=======================调用算法模板进行种群进化=============================="""
NDSet = myAlgorithm.run() # 执行算法模板，得到非支配种群
NDSet.save()              # 把结果保存到文件中
# 输出
print('用时：%f 秒'%(myAlgorithm.passTime))
print('评价次数：%d 次'%(myAlgorithm.evalsNum))
print('非支配个体数：%d 个'%(NDSet.sizes))
print('单位时间找到帕累托前沿点个数：%d 个'%(int(NDSet.sizes // myAlgorithm.passTime)))

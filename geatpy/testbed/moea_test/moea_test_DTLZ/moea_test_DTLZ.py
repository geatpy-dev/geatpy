# -*- coding: utf-8 -*-
import geatpy as ea # import geatpy

"""==================================实例化问题对象================================"""
problemName = 'DTLZ1' # 目标函数名
M = 3 # 设置目标维数
fileName = problemName # 这里因为目标函数写在与之同名的文件里，所以文件名也是目标函数名
MyProblem = getattr(__import__(fileName), problemName) # 获得自定义问题类
problem = MyProblem(M) # 生成问题对象

"""==================================种群设置================================"""
Encoding = 'R'             # 编码方式
conordis = 0               # 表示染色体解码后得到的变量是连续的
NIND = 100                 # 种群规模
Field = ea.crtfld(Encoding, conordis, problem.ranges, problem.borders) # 创建区域描述器
population = ea.Population(Encoding, conordis, Field, NIND) # 实例化种群对象（此时种群还没被真正初始化）

"""==================================算法参数设置================================"""
myAlgorithm = ea.moea_NSGA3_templet(problem, population) # 实例化一个算法模板对象
myAlgorithm.MAXGEN = 500 # 最大遗传代数
"""=======================调用算法模板进行种群进化=============================="""
NDSet = myAlgorithm.run() # 执行算法模板，得到帕累托最优解集NDSet
# 输出
print('用时：%s 秒'%(myAlgorithm.passTime))
print('非支配个体数：%s 个'%(NDSet.sizes))
print('单位时间找到帕累托前沿点个数：%s 个'%(int(NDSet.sizes // myAlgorithm.passTime)))
# 计算指标
PF = problem.calBest() # 计算真实前沿
if PF is not None:
    GD = ea.indicator.GD(NDSet.ObjV, PF) # 计算GD指标
    IGD = ea.indicator.IGD(NDSet.ObjV, PF) # 计算IGD指标
    HV = ea.indicator.HV(NDSet.ObjV, PF)
    Space = ea.indicator.spacing(NDSet.ObjV)
    print('GD',GD)
    print('IGD',IGD)
    print('HV', HV)
    print('Space', Space)

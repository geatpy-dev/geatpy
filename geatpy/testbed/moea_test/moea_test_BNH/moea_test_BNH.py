# -*- coding: utf-8 -*-
import geatpy as ea # import geatpy

"""================================实例化问题对象============================="""
problemName = 'BNH'         # 问题名称
fileName = problemName      # 这里因为目标函数写在与之同名的文件里，所以文件名也是问题名称
MyProblem = getattr(__import__(fileName), problemName) # 获得自定义问题类
problem = MyProblem()       # 生成问题对象
"""==================================种群设置================================"""
Encoding = 'RI'             # 编码方式
NIND = 200                  # 种群规模
Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders) # 创建区域描述器
population = ea.Population(Encoding, Field, NIND) # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
"""================================算法参数设置==============================="""
myAlgorithm = ea.moea_NSGA2_templet(problem, population) # 实例化一个算法模板对象
myAlgorithm.MAXGEN = 50     # 最大进化代数（即从初始种群开始最多进化多少代）
myAlgorithm.drawing = 1     # 设置绘图方式（0：不绘图；1：绘制结果图；2：绘制过程动画）
"""===========================调用算法模板进行种群进化===========================
调用run执行算法模板，得到帕累托最优解集NDSet。NDSet是一个种群类Population的对象。
NDSet.ObjV为最优解个体的目标函数值；NDSet.Phen为对应的决策变量值。
详见Population.py中关于种群类的定义。
"""
NDSet = myAlgorithm.run()   # 执行算法模板，得到非支配种群
NDSet.save()                # 把结果保存到文件中
# 输出
print('用时：%f 秒'%(myAlgorithm.passTime))
print('评价次数：%d 次'%(myAlgorithm.evalsNum))
print('非支配个体数：%d 个'%(NDSet.sizes))
print('单位时间找到帕累托前沿点个数：%d 个'%(int(NDSet.sizes // myAlgorithm.passTime)))
# 计算指标
PF = problem.getBest()      # 获取真实前沿，详见Problem.py中关于Problem类的定义
if PF is not None and NDSet.sizes != 0:
    GD = ea.indicator.GD(NDSet.ObjV, PF)       # 计算GD指标
    IGD = ea.indicator.IGD(NDSet.ObjV, PF)     # 计算IGD指标
    HV = ea.indicator.HV(NDSet.ObjV, PF)       # 计算HV指标
    Spacing = ea.indicator.Spacing(NDSet.ObjV) # 计算Spacing指标
    print('GD',GD)
    print('IGD',IGD)
    print('HV', HV)
    print('Spacing', Spacing)
"""=============================进化过程指标追踪分析============================"""
if PF is not None:
    metricName = [['IGD'], ['HV']]
    [NDSet_trace, Metrics] = ea.indicator.moea_tracking(myAlgorithm.pop_trace, PF, metricName, problem.maxormins)
    # 绘制指标追踪分析图
    ea.trcplot(Metrics, labels = metricName, titles = metricName)

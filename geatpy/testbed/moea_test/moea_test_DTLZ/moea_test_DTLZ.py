# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea  # import geatpy

if __name__ == '__main__':
    """==============================实例化问题对象============================="""
    problemName = 'DTLZ1'  # 问题名称
    M = 3  # 设置目标维数
    fileName = problemName  # 这里因为目标函数写在与之同名的文件里，所以文件名也是问题名称
    MyProblem = getattr(__import__(fileName), problemName)  # 获得自定义问题类
    problem = MyProblem(M)  # 生成问题对象
    """=================================种群设置==============================="""
    Encoding = 'RI'  # 编码方式
    NIND = 100  # 种群规模
    Field = ea.crtfld(Encoding, list(problem.varTypes), problem.ranges, problem.borders)  # 创建区域描述器
    population = ea.Population(Encoding, Field, NIND)  # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
    """===============================算法参数设置=============================="""
    myAlgorithm = ea.moea_NSGA3_templet(problem, population)  # 实例化一个算法模板对象
    myAlgorithm.MAXGEN = 500  # 最大进化代数
    myAlgorithm.logTras = 1  # 设置每多少代记录日志，若设置成0则表示不记录日志
    myAlgorithm.verbose = True  # 设置是否打印输出日志信息
    myAlgorithm.drawing = 1  # 设置绘图方式（0：不绘图；1：绘制结果图；2：绘制目标空间过程动画；3：绘制决策空间过程动画）
    """==========================调用算法模板进行种群进化=========================
    调用run执行算法模板，得到帕累托最优解集NDSet以及最后一代种群。NDSet是一个种群类Population的对象。
    NDSet.ObjV为最优解个体的目标函数值；NDSet.Phen为对应的决策变量值。
    详见Population.py中关于种群类的定义。
    """
    [NDSet, population] = myAlgorithm.run()  # 执行算法模板，得到非支配种群以及最后一代种群
    NDSet.save()  # 把非支配种群的信息保存到文件中
    """==================================输出结果=============================="""
    print('用时：%f 秒' % myAlgorithm.passTime)
    print('评价次数：%d 次' % myAlgorithm.evalsNum)
    print('非支配个体数：%d 个' % NDSet.sizes) if NDSet.sizes != 0 else print('没有找到可行解！')
    if myAlgorithm.log is not None and NDSet.sizes != 0:
        print('GD', myAlgorithm.log['gd'][-1])
        print('IGD', myAlgorithm.log['igd'][-1])
        print('HV', myAlgorithm.log['hv'][-1])
        print('Spacing', myAlgorithm.log['spacing'][-1])
        """=========================进化过程指标追踪分析========================="""
        metricName = [['igd'], ['hv']]
        Metrics = np.array([myAlgorithm.log[metricName[i][0]] for i in range(len(metricName))]).T
        # 绘制指标追踪分析图
        ea.trcplot(Metrics, labels=metricName, titles=metricName)

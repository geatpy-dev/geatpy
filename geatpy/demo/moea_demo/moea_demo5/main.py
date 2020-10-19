# -*- coding: utf-8 -*-
import geatpy as ea  # import geatpy
from MyProblem import MyProblem  # 导入自定义问题接口

if __name__ == '__main__':
    """===============================实例化问题对象============================"""
    problem = MyProblem()  # 生成问题对象
    """==================================种群设置==============================="""
    Encoding = 'RI'  # 编码方式
    NIND = 50  # 种群规模
    Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders)  # 创建区域描述器
    population = ea.Population(Encoding, Field, NIND)  # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
    """================================算法参数设置============================="""
    myAlgorithm = ea.moea_awGA_templet(problem, population)  # 实例化一个算法模板对象
    myAlgorithm.MAXGEN = 20  # 最大进化代数
    myAlgorithm.logTras = 0  # 设置每多少代记录日志，若设置成0则表示不记录日志
    myAlgorithm.verbose = False  # 设置是否打印输出日志信息
    myAlgorithm.drawing = 0  # 设置绘图方式（0：不绘图；1：绘制结果图；2：绘制目标空间过程动画；3：绘制决策空间过程动画）
    """==========================调用算法模板进行种群进化========================"""
    prophetPop, pop = myAlgorithm.run()  # 执行算法模板，得到帕累托最优解集，作为一个“先知种群”
    """============================插入先验知识进行进化=========================="""
    myAlgorithm = ea.moea_NSGA2_templet(problem, population)  # 实例化一个算法模板对象
    myAlgorithm.MAXGEN = 50  # 最大进化代数
    myAlgorithm.logTras = 1
    myAlgorithm.verbose = True
    myAlgorithm.drawing = 1
    [NDSet, population] = myAlgorithm.run(prophetPop)  # 插入“先知种群”，执行算法模板
    NDSet.save()  # 把非支配种群的信息保存到文件中
    """==================================输出结果=============================="""
    print('用时：%s 秒' % myAlgorithm.passTime)
    print('非支配个体数：%d 个' % NDSet.sizes) if NDSet.sizes != 0 else print('没有找到可行解！')

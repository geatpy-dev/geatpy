# -*- coding: utf-8 -*-
import geatpy as ea  # import geatpy
from MyProblem import MyProblem  # 导入自定义问题接口

if __name__ == '__main__':
    """===============================实例化问题对象============================"""
    problem = MyProblem()  # 生成问题对象
    """==================================种群设置=============================="""
    NIND = 30  # 种群规模
    # 创建区域描述器，这里需要创建两个，前1个变量用RI编码，剩余变量用排列编码
    Encodings = ['RI', 'P']  # 编码方式列表
    Field1 = ea.crtfld(Encodings[0], problem.varTypes[:1], problem.ranges[:, :1], problem.borders[:, :1])  # 创建区域描述器
    Field2 = ea.crtfld(Encodings[1], problem.varTypes[1:], problem.ranges[:, 1:], problem.borders[:, 1:])
    Fields = [Field1, Field2]
    population = ea.PsyPopulation(Encodings, Fields, NIND)  # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
    """================================算法参数设置============================"""
    myAlgorithm = ea.moea_psy_NSGA2_templet(problem, population)  # 实例化一个算法模板对象
    myAlgorithm.MAXGEN = 200  # 最大进化代数
    myAlgorithm.logTras = 0  # 设置每多少代记录日志，若设置成0则表示不记录日志
    myAlgorithm.verbose = False  # 设置是否打印输出日志信息
    myAlgorithm.drawing = 1  # 设置绘图方式（0：不绘图；1：绘制结果图；2：绘制目标空间过程动画；3：绘制决策空间过程动画）
    """==========================调用算法模板进行种群进化=========================
    调用run执行算法模板，得到帕累托最优解集NDSet以及最后一代种群。NDSet是一个种群类Population的对象。
    NDSet.ObjV为最优解个体的目标函数值；NDSet.Phen为对应的决策变量值。
    详见Population.py中关于种群类的定义。
    """
    [NDSet, population] = myAlgorithm.run()  # 执行算法模板，得到非支配种群以及最后一代种群
    NDSet.save()  # 把非支配种群的信息保存到文件中
    """==================================输出结果=============================="""
    print('用时：%s 秒' % myAlgorithm.passTime)
    print('非支配个体数：%d 个' % NDSet.sizes) if NDSet.sizes != 0 else print('没有找到可行解！')

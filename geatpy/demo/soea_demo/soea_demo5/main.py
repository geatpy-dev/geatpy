# -*- coding: utf-8 -*-
import geatpy as ea  # import geatpy
from MyProblem import MyProblem  # 导入自定义问题接口

if __name__ == '__main__':
    """================================实例化问题对象==========================="""
    problem = MyProblem()  # 生成问题对象
    """==================================种群设置=============================="""
    NIND = 40  # 种群规模
    # 创建区域描述器，这里需要创建两个，前2个变量用RI编码，剩余变量用排列编码
    Encodings = ['RI', 'P']
    Field1 = ea.crtfld(Encodings[0], problem.varTypes[:2], problem.ranges[:, :2], problem.borders[:, :2])
    Field2 = ea.crtfld(Encodings[1], problem.varTypes[2:], problem.ranges[:, 2:], problem.borders[:, 2:])
    Fields = [Field1, Field2]
    population = ea.PsyPopulation(Encodings, Fields, NIND)  # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
    """================================算法参数设置============================="""
    myAlgorithm = ea.soea_psy_EGA_templet(problem, population)  # 实例化一个算法模板对象
    myAlgorithm.MAXGEN = 25  # 最大进化代数
    myAlgorithm.logTras = 1  # 设置每隔多少代记录日志，若设置成0则表示不记录日志
    myAlgorithm.verbose = True  # 设置是否打印输出日志信息
    myAlgorithm.drawing = 1  # 设置绘图方式（0：不绘图；1：绘制结果图；2：绘制目标空间过程动画；3：绘制决策空间过程动画）
    """===========================调用算法模板进行种群进化========================"""
    [BestIndi, population] = myAlgorithm.run()  # 执行算法模板，得到最优个体以及最后一代种群
    BestIndi.save()  # 把最优个体的信息保存到文件中
    """==================================输出结果=============================="""
    print('评价次数：%s' % myAlgorithm.evalsNum)
    print('时间已过 %s 秒' % myAlgorithm.passTime)
    if BestIndi.sizes != 0:
        print('最优的目标函数值为：%s' % (BestIndi.ObjV[0][0]))
        print('最优的控制变量值为：')
        for i in range(BestIndi.Phen.shape[1]):
            print(BestIndi.Phen[0, i])
    else:
        print('没找到可行解。')

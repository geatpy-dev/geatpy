# -*- coding: utf-8 -*-
import geatpy as ea  # import geatpy

if __name__ == '__main__':
    """===============================实例化问题对象============================"""
    problemName = 'Beale'  # 问题名称
    fileName = problemName  # 这里因为目标函数写在与之同名的文件里，所以文件名也是问题名称
    MyProblem = getattr(__import__(fileName), problemName)  # 获得自定义问题类
    problem = MyProblem()  # 生成问题对象
    """==================================种群设置=============================="""
    Encoding = 'RI'  # 编码方式
    NIND = 20  # 种群规模
    Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders)  # 创建区域描述器
    population = ea.Population(Encoding, Field, NIND)  # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
    """================================算法参数设置============================="""
    myAlgorithm = ea.soea_DE_rand_1_bin_templet(problem, population)  # 实例化一个算法模板对象
    myAlgorithm.MAXGEN = 100  # 最大进化代数
    myAlgorithm.mutOper.F = 0.5  # 差分进化中的参数F
    myAlgorithm.recOper.XOVR = 0.2  # 重组概率
    myAlgorithm.logTras = 1  # 设置每隔多少代记录日志，若设置成0则表示不记录日志
    myAlgorithm.verbose = True  # 设置是否打印输出日志信息
    myAlgorithm.drawing = 1  # 设置绘图方式（0：不绘图；1：绘制结果图；2：绘制目标空间过程动画；3：绘制决策空间过程动画）
    """===========================调用算法模板进行种群进化======================="""
    [BestIndi, population] = myAlgorithm.run()  # 执行算法模板，得到最优个体以及最后一代种群
    BestIndi.save()  # 把最优个体的信息保存到文件中
    """==================================输出结果=============================="""
    print('用时：%f 秒' % myAlgorithm.passTime)
    print('评价次数：%d 次' % myAlgorithm.evalsNum)
    if BestIndi.sizes != 0:
        print('最优的目标函数值为：%s' % BestIndi.ObjV[0][0])
        print('最优的控制变量值为：')
        for i in range(BestIndi.Phen.shape[1]):
            print(BestIndi.Phen[0, i])
    else:
        print('没找到可行解。')

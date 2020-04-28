# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea   # import geatpy
from MyProblem import MyProblem # 导入自定义问题接口

if __name__ == '__main__':
    """================================实例化问题对象==========================="""
    problem = MyProblem() # 生成问题对象
    """==================================种群设置==============================="""
    Encoding = 'RI'       # 编码方式
    NIND = 20             # 种群规模
    Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders) # 创建区域描述器
    population = ea.Population(Encoding, Field, NIND) # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
    """================================算法参数设置============================="""
    myAlgorithm = ea.soea_DE_currentToBest_1_bin_templet(problem, population) # 实例化一个算法模板对象
    myAlgorithm.MAXGEN = 400 # 最大进化代数
    myAlgorithm.mutOper.F = 0.7 # 差分进化中的参数F
    myAlgorithm.recOper.XOVR = 0.7 # 重组概率
    myAlgorithm.drawing = 1
    """===========================根据先验知识创建先知种群======================="""
    prophetChrom = np.array([[0.4, 0.2, 0.4]]) # 假设已知[0.4, 0.2, 0.4]为一条比较优秀的染色体
    prophetPop = ea.Population(Encoding, Field, 1, prophetChrom) # 实例化种群对象（设置个体数为1）
    myAlgorithm.call_aimFunc(prophetPop) # 计算先知种群的目标函数值及约束（假如有约束）
    """===========================调用算法模板进行种群进化======================="""
    [population, obj_trace, var_trace] = myAlgorithm.run(prophetPop) # 执行算法模板，插入先知种群
    population.save() # 把最后一代种群的信息保存到文件中
    # 输出结果
    best_gen = np.argmin(problem.maxormins * obj_trace[:, 1]) # 记录最优种群个体是在哪一代
    best_ObjV = obj_trace[best_gen, 1]
    print('最优的目标函数值为：%s'%(best_ObjV))
    print('最优的决策变量值为：')
    for i in range(var_trace.shape[1]):
        print(var_trace[best_gen, i])
    print('有效进化代数：%s'%(obj_trace.shape[0]))
    print('最优的一代是第 %s 代'%(best_gen + 1))
    print('评价次数：%s'%(myAlgorithm.evalsNum))
    print('时间已过 %s 秒'%(myAlgorithm.passTime))

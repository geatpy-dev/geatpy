# -*- coding: utf-8 -*-

import geatpy as ea
import numpy as np

"""

本案例展示了如何快速创建问题对象、快速调用算法求解一个含两个约束的单目标优化问题。
与soea_quick_start_evalVars案例不同的是，本案例通过装饰器single标记目标函数evalVars，使得目标函数只传入一组决策变量。
此时目标函数只需计算这组变量对应的目标函数值以及违反约束程度值。
用法：
    ObjV = evalVars(Vars) -- 若无约束
    或
    ObjV, CV = evalVars(Vars) -- 若有约束
注意：
    加了装饰器single后，evalVars的传入参数Vars不再是Numpy ndarray二维数组，而是Numpy ndarray一维数组。
    在设置返回值时：
        返回值ObjV可以赋值为一个Numpy ndarray一维数组，也可以是一个标量（当只有一个优化目标时）。
        返回值CV可以赋值为一个Numpy ndarray一维数组，也可以是一个标量（当只有一个优化目标时）。
    值得注意的是：
        通过调用evalVars(Vars)得到的返回值会被single装饰器修正为Numpy ndarray二维数组。
        即ObjV = evalVars(Vars)得到的ObjV或：ObjV, CV = evalVars(Vars)得到的ObjV和CV，均为Numpy ndarray二维数组。

"""

if __name__ == '__main__':
    # 构建问题
    r = 1  # 目标函数需要用到的额外数据
    @ea.Problem.single
    def evalVars(Vars):  # 定义目标函数（含约束）
        f = np.sum((Vars - r) ** 2)  # 计算目标函数值
        x1 = Vars[0]
        x2 = Vars[1]
        CV = np.array([(x1 - 0.5)**2 - 0.25,
                       (x2 - 1)**2 - 1])  # 计算违反约束程度
        return f, CV

    problem = ea.Problem(name='soea quick start demo',
                         M=1,  # 目标维数
                         maxormins=[1],  # 目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标
                         Dim=5,  # 决策变量维数
                         varTypes=[0, 0, 1, 1, 1],  # 决策变量的类型列表，0：实数；1：整数
                         lb=[-1, 1, 2, 1, 0],  # 决策变量下界
                         ub=[1, 4, 5, 2, 1],  # 决策变量上界
                         evalVars=evalVars)
    # 构建算法
    algorithm = ea.soea_SEGA_templet(problem,
                                     ea.Population(Encoding='RI', NIND=20),
                                     MAXGEN=50,  # 最大进化代数。
                                     logTras=1,  # 表示每隔多少代记录一次日志信息，0表示不记录。
                                     trappedValue=1e-6,  # 单目标优化陷入停滞的判断阈值。
                                     maxTrappedCount=10)  # 进化停滞计数器最大上限值。
    # 求解
    res = ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=True, drawLog=False, saveFlag=True)
    print(res)

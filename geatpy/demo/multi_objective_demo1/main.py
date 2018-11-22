# -*- coding: utf-8 -*-
"""
执行脚本main.py
描述：
    该demo是展示如何计算带约束的两个最小化目标的帕累托前沿
    两个目标为：
    min -25 * (x1 - 2)**2 - (x2 - 2)**2 - (x3 - 1)**2 - (x4 - 4)**2 - (x5 - 1)**2
    min (x1 - 1)**2 + (x2 - 1)**2 + (x3 - 1)**2 + (x4 - 1)**2 + (x5 - 1)**2
    约束条件为：
    x1 + x2 >= 2
    x1 + x2 <= 6
    x1 - x2 >= -2
    x1 - 3*x2 <= 2
    4 - (x3 - 3)**2 - x4 >= 0
    (x5 - 3)**2 + x4 - 4 >= 0
    其中目标函数和约束条件写在aimfuc.py文件中，适应度罚函数写在罚函数文件punishing.py中
"""

import geatpy as ga # 导入geatpy库
from multimin import multimin # 导入自定义的进化算法模板

# 获取函数接口地址
AIM_M = __import__('aimfuc')
PUN_M = __import__('punishing')
"""============================变量设置============================"""
NVAR = 50 # 变量个数
Base = 11 # 变量取值
"""========================遗传算法参数设置========================="""
NIND = 50;               # 种群规模
MAXGEN = 200             # 最大遗传代数
GGAP = 0.8               # 代沟：子代与父代的重复率为(1-GGAP)
selectStyle = 'rws'      # 遗传算法的选择方式设为"rws"——轮盘赌选择
recombinStyle = 'xovdp'  # 遗传算法的重组方式，设为两点交叉
recopt = 0.9             # 交叉概率
pm = 0.1                 # 变异概率
SUBPOP = 1               # 设置种群数为1
maxormin = 1             # 设置标记表明这是最小化目标
"""=======================调用编程模板进行种群进化==================="""
# 调用编程模板进行种群进化，得到种群进化和变量的追踪器以及运行时间
[ObjV, NDSet, NDSetObjV, times] = multimin(AIM_M, 'aimfuc', PUN_M, 'punishing', NIND, NVAR, Base, MAXGEN, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, maxormin)
"""=========================绘图及输出结果========================="""
ga.frontplot(NDSetObjV, True)
print('用时：', times, '秒')
print('帕累托最优解个数：', NDSet.shape[0])
print('平均每秒找到的帕累托前沿个数: ', int(NDSet.shape[0] / times))

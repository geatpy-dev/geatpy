"""最小化目标函数值问题求解程序执行脚本main.py"""
import numpy as np
import geatpy as ga # 导入geatpy库
from multimin import multimin # 导入自定义的编程模板
from aimfuc import aimfuc # 导入自定义的目标函数接口

# 获取函数接口地址
AIM_M = __import__('aimfuc')

"""============================变量设置============================"""
NVAR = 50 # 变量个数
Base = 11 # 变量取值
"""========================遗传算法参数设置========================="""
NIND = 50;               # 种群规模
MAXGEN = 2000;             # 最大遗传代数
GGAP = 0.8;              # 代沟：子代与父代的重复率为(1-GGAP)
selectStyle = 'rws'      # 遗传算法的选择方式设为"rws"——轮盘赌选择
recombinStyle = 'xovdp'  # 遗传算法的重组方式，设为两点交叉
recopt = 0.9             # 交叉概率
pm = 0.1                # 变异概率
SUBPOP = 1               # 设置种群数为1
maxormin = 1             # 设置标记表明这是最小化目标
"""=======================调用编程模板进行种群进化==================="""
# 调用编程模板进行种群进化，得到种群进化和变量的追踪器以及运行时间
[ObjV, NDSet, times] = multimin(AIM_M, 'aimfuc', NIND, NVAR, Base, MAXGEN, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, maxormin)
"""=========================绘图及输出结果========================="""
ga.frontplot(NDSet, True)
print('用时：', times, '秒')
print(NDSet.shape)
print('平均每秒找到的帕累托前沿个数: ', int(NDSet.shape[0] / times))
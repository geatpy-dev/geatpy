"""
执行脚本main.py
描述：
    该demo是展示如何求解二元二次方程的整数解：
    x + y = xy
    本案例通过自定义算法模板"mintemp1.py"来解决该问题，
    mintemp1这个算法模板可以搜索出具有相同目标函数值的不同的解
    先预估解在[-100, 100]之间，再逐渐缩小范围，最后在[-5,5]之间进行搜索
    其中目标函数和约束条件写在aimfuc.py文件中
"""
import numpy as np
import geatpy as ga # 导入geatpy库
from mintemp1 import mintemp1 # 导入自定义的编程模板

# 获取函数接口地址
AIM_M = __import__('aimfuc')        # 目标函数
"""============================变量设置============================"""
# 变量设置
ranges = np.array([[-5, -5],[5, 5]]) # 生成自变量的范围矩阵
borders = np.ones((2, 2)) # 生成自变量的边界矩阵
FieldDR = ga.crtfld(ranges, borders) # 生成区域描述器
"""=======================调用编程模板进行种群进化==================="""
[pop_trace, solutions, times] = mintemp1(AIM_M, 'aimfuc', None, None, FieldDR, problem = 'I', maxormin = 1, MAXGEN = 50, NIND = 50, SUBPOP = 1, GGAP = 0.9, selectStyle = 'sus', recombinStyle = 'xovdprs', recopt = 0.9, pm = 0.5, distribute = True, drawing = 1)

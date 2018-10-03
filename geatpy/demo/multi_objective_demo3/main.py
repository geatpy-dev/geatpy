# -*- coding: utf-8 -*-
"""
执行脚本main.py
描述：
    该demo是展示如何计算带约束的两个最大化目标的帕累托前沿
    本案例调用q_sorted_new_templet算法模板解决该问题
    其中为了方便展示，把目标函数和罚函数均和执行脚本写在同一个文件里，建议分开写在不同的文件中。
"""

import numpy as np
import geatpy as ga

# 注意：不建议把目标函数放在执行脚本内，建议放在另一个文件中
def aimfuc(x): # 定义目标函数
    x1 = x[:, [0]]
    fun1 = -x1**2
    fun2 = -(x1 - 2)**2
    return np.hstack([fun1, fun2]) # 对矩阵进行转置使得目标函数矩阵符合Geatpy数据结构

def punishing(x, FitnV): # 定义罚函数
    x1 = x[:, [0]]
    idx1 = np.where(x1**2 - 2.5 * x1 + 1.5 < 0)[0]
    FitnV[idx1] = 0
    return [FitnV, idx1] # 返回新的适应度以及非可行解的下标矩阵

if __name__ == "__main__":
    AIM_M = __import__('main') # 获取函数所在文件的地址
    PUN_M = __import__('main') # 获取罚函数所在文件的地址
    # 变量设置
    ranges = np.array([[-10], [10]])  # 生成自变量的范围矩阵
    borders = np.array([[1], [1]])   # 生成自变量的边界矩阵（1表示变量的区间是闭区间）
    precisions = [1] # 因为变量范围都是闭区间，而且nsga2_templet编程模板采用的是实数编码，因此precisions不起作用，但要设置成大于0的任意值
    FieldDR = ga.crtfld(ranges, borders, precisions) # 生成区域描述器
    # 调用编程模板
    [ObjV, NDSet, NDSetObjV, times] = ga.q_sorted_new_templet(AIM_M, 'aimfuc', PUN_M, 'punishing', FieldDR, 'R', maxormin = -1, MAXGEN = 1000, MAXSIZE = 500, NIND = 50, SUBPOP = 1, GGAP = 1, selectStyle = 'tour', recombinStyle = 'xovdp', recopt = 0.9, pm = 0.6, distribute = True, drawing = 1)
    
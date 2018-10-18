# -*- coding: utf-8 -*-
"""
执行脚本main.py
描述：
    该demo是展示如何计算无约束条件的两个最小化目标的帕累托前沿。
    其中为了方便，展示如何直接将目标函数和执行脚本写在同一个文件里，
    类似地，适应度罚函数也可以和执行脚本写在同一个文件里，但最好是分开写在另一个文件中。
    本案例调用了“moea_nsga2_templet”算法模板，其详细用法可利用help命令查看，或是在github下载并查看源码。
    调用算法模板时可以设置drawing=2，此时算法模板将在种群进化过程中绘制动画，但注意执行前要在Python控制台执行命令matplotlib qt5。
"""
import numpy as np
import geatpy as ga

# 注意：不建议把目标函数放在执行脚本内，建议放在另一个文件中
def aimfuc(x,LegV): # 定义目标函数
    x1 = x[:, 0]; x2 = x[:, 1]
    fun1 = x1**4-10*x1**2+x1*x2+x2**4-x1**2*x2**2
    fun2 = x2**4-x1**2*x2**2+x1**4+x1*x2
    return [np.vstack([fun1, fun2]).T ,LegV] # 对矩阵进行转置使得目标函数矩阵符合Geatpy数据结构

if __name__ == "__main__":
    AIM_M = __import__('main') # 获取函数接口所在文件的地址
    # 变量设置
    ranges = np.array([[-5, -5], [5, 5]])  # 生成自变量的范围矩阵
    borders = np.array([[1, 1], [1, 1]])   # 生成自变量的边界矩阵（1表示变量的区间是闭区间）
    precisions = [1, 1] # 根据crtfld的函数特性，这里需要设置精度为任意正值，否则在生成区域描述器时会默认为整数编码，并对变量范围作出一定调整
    FieldDR = ga.crtfld(ranges, borders, precisions) # 生成区域描述器
    # 调用编程模板
    [ObjV, NDSet, NDSetObjV, times] = ga.moea_nsga2_templet(AIM_M, 'aimfuc', None, None, FieldDR, 'R', maxormin = 1, MAXGEN = 500, MAXSIZE = 200, NIND = 100, SUBPOP = 1, GGAP = 1, selectStyle = 'tour', recombinStyle = 'xovdp', recopt = 0.9, pm = 0.6, distribute = True, drawing = 1)
    
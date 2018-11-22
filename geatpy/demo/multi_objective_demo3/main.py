# -*- coding: utf-8 -*-
"""
执行脚本main.py
描述：
    该demo是展示如何计算带约束的两个最大化目标的帕累托前沿
    min x1**2
    min (x1 - 2)**2
    s.t. 
    x1 ∈ [-10,10]
    x1**2 - 2.5 * x1 + 1.5 >= 0
    其中为了方便展示，把目标函数aimfuc和适应度罚函数punishing均和执行脚本写在同一个文件里，建议分开写在不同的文件中。
    本案例调用了“moea_q_sorted_new_templet”算法模板，其详细用法可利用help命令查看，或是在github下载并查看源码。
    调用算法模板时可以设置drawing=2，此时算法模板将在种群进化过程中绘制动画，但注意执行前要在Python控制台执行命令matplotlib qt5。
"""

import numpy as np
import geatpy as ga

# 注意：不建议把目标函数放在执行脚本内，建议放在另一个文件中
def aimfuc(x, LegV): # 定义目标函数
    x1 = x[:, [0]]
    fun1 = x1**2
    fun2 = (x1 - 2)**2
    # 约束条件
    exIdx = np.where(x1**2 - 2.5 * x1 + 1.5 < 0)[0] # 获取不满足约束条件的个体在种群中的下标
    # 惩罚非可行解方法1：修改非可行个体的目标函数值的方法，对非可行解对应的目标函数值作出惩罚，而不标记其为非可行解
    # 不标记非可行解的做法比较危险，假如处理不慎，可能会导致以下错误情况：
    # 某一代的种群全是非可行解而其组成的“不合法的帕累托最优解集”恰好支配了其他所有的可行解。
    # 此时进化记录器会受到欺骗，而得出一个实际上是非可行解的搜索结果。
    # 一个比较好的处理方法是：要让其绝对地比其他任何可行解的目标函数值都要大或小（看是最小化目标还是最大化目标）。
#    fun1[exIdx] = np.max(fun1) # 此处修改成max(ObjV1)是危险的，因为无法保证修改后能被其他可行解支配
#    fun2[exIdx] = np.max(fun2)
    
    # 惩罚非可行解的方法2：修改非可行个体的适应度的方法，配合punishing罚函数接口，先标记非可行解，然后在punishing中对其适应度FitnV作出惩罚
    # (也可以不写punishing，因为算法模板中在一些与适应度有关的计算时已传入LegV，此时内核函数会自动处理非可行解的适应度)
    # (写punishing实质上是进一步惩罚非可行解的适应度值)
    LegV[exIdx] = 0 # 标记非可行解对应的可行性列向量中元素的值为0
    return [np.hstack([fun1, fun2]), LegV] # 对矩阵进行转置使得目标函数矩阵符合Geatpy数据结构

def punishing(LegV, FitnV): # 定义罚函数
    FitnV[np.where(LegV == 0)[0]] = 0
    return FitnV # 返回新的适应度

if __name__ == "__main__":
    AIM_M = __import__('main') # 获取函数所在文件的地址
    PUN_M = __import__('main') # 获取罚函数所在文件的地址
    # 变量设置
    ranges = np.array([[-10], [10]])  # 生成自变量的范围矩阵
    borders = np.array([[1], [1]])   # 生成自变量的边界矩阵（1表示变量的区间是闭区间）
    precisions = [1] # 根据crtfld的函数特性，这里需要设置精度为任意正值，否则在生成区域描述器时会默认为整数编码，并对变量范围作出一定调整
    FieldDR = ga.crtfld(ranges, borders, precisions) # 生成区域描述器
    # 调用编程模板
    [ObjV, NDSet, NDSetObjV, times] = ga.moea_q_sorted_new_templet(AIM_M, 'aimfuc', PUN_M, 'punishing', FieldDR, 'R', maxormin = 1, MAXGEN = 100, MAXSIZE = 500, NIND = 50, SUBPOP = 1, GGAP = 1, selectStyle = 'etour', recombinStyle = 'xovdprs', recopt = 0.9, pm = 0.1, distribute = True, drawing = 1)
    
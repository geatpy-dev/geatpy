# -*- coding: utf-8 -*-
"""
执行脚本main.py
描述：
    该demo是展示如何计算带等式约束的单目标优化问题:
        max 4 * x1 + 2 * x2 + x3
        s.t.
            2 * x1 + x2 <= 1
            x1 + 2 * x3 > 2
            x1 + x2 + x3 == 1
            0 <= x1 <= 1
            0 <= x2 <= 1
            0 <= x3 <= 2
    其中目标函数和约束条件写在aimfuc.py文件中，适应度罚函数写在罚函数文件punishing.py中
    本案例通过降维的方法，将等式约束化成了不等式约束，大大拓宽了可行解的空间，方便遗传算法求解
    此外，本案例展示了利用多种群竞争的进化算法模板sga_mpc_real_templet了解决该问题。
"""

import numpy as np
import geatpy as ga

# 获取函数接口地址
AIM_M = __import__('aimfuc')
PUN_M = __import__('punishing')
# 变量设置
x1 = [0, 1] # 自变量1的范围
x2 = [0, 1] # 自变量2的范围
b1 = [1, 1] # 自变量1是否包含下界
b2 = [1, 1] # 自变量2是否包含上界
ranges=np.vstack([x1, x2]).T # 生成自变量的范围矩阵
borders = np.vstack([b1, b2]).T # 生成自变量的边界矩阵
precisions = [2] * 2 # 在二进制/格雷码编码中代表自变量的编码精度，当控制变量是连续型时，根据crtfld参考资料，该变量只表示边界精度，故设置为一定的正数即可
newRanges = ga.meshrng(ranges, gridnum = 2) # 对控制变量范围进行网格化，网格边长为2
# 生成网格化后的区域描述器集合
FieldDRs = []
for i in range(len(newRanges)):
    FieldDRs.append(ga.crtfld(newRanges[i], borders, precisions))
# 调用编程模板
[pop_trace, var_trace, times] = ga.sga_mpc_real_templet(AIM_M, 'aimfuc', PUN_M,\
 'punishing', FieldDRs, problem = 'R', maxormin = -1, MAXGEN = 50, NIND = 50,\
 SUBPOP = 1, GGAP = 0.9, selectStyle = 'tour', recombinStyle = 'xovdprs',\
 recopt = 0.9, pm = 0.3, distribute = True, drawing = 1)

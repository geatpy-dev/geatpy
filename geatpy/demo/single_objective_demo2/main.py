# -*- coding: utf-8 -*-
"""
执行脚本main.py
描述：
    该demo是展示如何计算无约束的单目标优化问题
    本案例通过调用sga_new_code_templet算法模板来解决该问题
    其中目标函数写在aimfuc.py文件中
    本案例调用了“sga_new_code_templet”这个算法模板，其详细用法可利用help命令查看，或是在github下载并查看源码
    调用算法模板时可以设置drawing=2，此时算法模板将在种群进化过程中绘制动画，但注意执行前要在Python控制台执行命令matplotlib qt5。
"""

import numpy as np
import geatpy as ga

# 获取函数接口地址
AIM_M = __import__('aimfuc')
# 变量设置
x1 = [-3, 12.1] # 自变量1的范围
x2 = [4.1, 5.8] # 自变量2的范围
b1 = [1, 1] # 自变量1是否包含下界
b2 = [1, 1] # 自变量2是否包含上界
codes = [0, 0] # 自变量的编码方式，0表示采用标准二进制编码
precisions = [4, 4] # 在二进制/格雷码编码中代表自变量的编码精度，当控制变量是连续型时，根据crtfld参考资料，该变量只表示边界精度，故设置为一定的正数即可
scales = [0, 0] # 是否采用对数刻度
ranges=np.vstack([x1, x2]).T # 生成自变量的范围矩阵
borders = np.vstack([b1, b2]).T # 生成自变量的边界矩阵
# 生成区域描述器
FieldD = ga.crtfld(ranges, borders, precisions, codes, scales)

# 调用编程模板
[pop_trace, var_trace, times] = ga.sga_new_code_templet(AIM_M, 'aimfuc', None, None, FieldD, problem = 'R', maxormin = -1, MAXGEN = 1000, NIND = 100, SUBPOP = 1, GGAP = 0.8, selectStyle = 'sus', recombinStyle = 'xovdp', recopt = None, pm = None, distribute = True, drawing = 1)

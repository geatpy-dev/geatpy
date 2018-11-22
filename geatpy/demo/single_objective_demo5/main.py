# -*- coding: utf-8 -*-
"""
执行脚本main.py
描述：
    该demo是展示如何计算带约束的任务指派问题
    其中目标函数和约束条件写在aimfuc.py文件中
问题如下：
    设有5个人，5个任务。
    已知这5个人每小时工作要求的工资分别是1,2,3,4,5元，
    而这5个任务分别耗时1,2,3,4,5小时。
    此外，已知工人1无法完成第2和第4个任务；工人3无法完成第1和第4个任务。
    现要求给每个人分配去完成不同的任务，要求老板一共支付工人的工钱数最少。
    因为问题需要用排列编码的染色体来解决，因此本案例调用了“sga_new_permut_templet”这个算法模板，其详细用法可利用help命令查看，或是在github下载并查看源码
    调用算法模板时可以设置drawing=2，此时算法模板将在种群进化过程中绘制动画，但注意执行前要在Python控制台执行命令matplotlib qt5。
"""

import geatpy as ga

# 获取函数接口地址
AIM_M = __import__('aimfuc')
# 参数设置
NVAR = 5                  # 排列编码的染色体长度
VarLen = 5                # 排列集合的大小，等于5表示排列集合为{1,2,3,4,5}
# 调用编程模板，其中recombinStyle要设置为'xovpm',对于排列编码问题，要采用特殊的xovpm(部分匹配交叉)的重组方式
[pop_trace, var_trace, times] = ga.sga_new_permut_templet(AIM_M, 'aimfuc', None, None, NVAR, VarLen, maxormin = 1, MAXGEN = 100, NIND = 10, SUBPOP = 1, GGAP = 0.9, selectStyle = 'etour', recombinStyle = 'xovpm', recopt = 0.9, pm = 0.1, distribute = True, drawing = 1)

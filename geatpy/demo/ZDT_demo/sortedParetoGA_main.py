"""sortedParetoGA_main.py"""
import numpy as np
import geatpy as ga # 导入geatpy库
# 导入自定义的编程模板
from awGA_tmplet import awGA
from i_awGA_tmplet import i_awGA
from q_sorted_tmplet import q_sorted
from nsga2_tmplet import nsga2
from aimfuc import aimfuc # 导入自定义的目标函数接口

# 获取函数接口地址
AIM_M = __import__('aimfuc')
AIM_F = 'aimfuc'

"""============================变量设置============================"""
ranges = np.vstack([np.zeros((1,30)), np.ones((1,30))])   # 生成自变量的范围矩阵
borders = np.vstack([np.ones((1,30)), np.ones((1,30))])       # 生成自变量的边界矩阵
precisions = [4] * 30              # 自变量的编码精度
"""========================遗传算法参数设置========================="""
NIND = 50              # 种群规模
MAXGEN = 500             # 最大遗传代数
GGAP = 1;              # 代沟：子代与父代的重复率为(1-GGAP)
selectStyle = 'tour'      # 遗传算法的选择方式
recombinStyle = 'xovdprs'  # 遗传算法的重组方式，设为两点交叉
recopt = 0.9             # 交叉概率
pm = None                   # 变异概率
SUBPOP = 1               # 设置种群数为1f
maxormin = 1             # 设置标记表明这是最小化目标
MAXSIZE = 1000                 # 帕累托最优集最大个数
"""=======================调用编程模板进行种群进化==================="""
[ObjV, NDSet, times] = ga.q_sorted_templet(AIM_M, AIM_F, None, None, ranges, borders, precisions, maxormin, MAXGEN, MAXSIZE, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, drawing = 1)


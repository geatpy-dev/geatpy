"""sortedParetoGA_main.py"""
import numpy as np
import geatpy as ga # import geatpy

# 获取函数接口地址
AIM_M = __import__('aimfuc')

"""============================变量设置============================"""
ranges = np.vstack([np.zeros((1,6)), np.ones((1,6))])   # 生成自变量的范围矩阵
borders = np.vstack([np.ones((1,6)), np.ones((1,6))])       # 生成自变量的边界矩阵
precisions = [4] * 30              # 自变量的编码精度
"""=======================调用编程模板进行种群进化==================="""
[ObjV, NDSet, times] = ga.nsga2_templet(AIM_M, 'aimfuc',None, None, ranges, borders, precisions, maxormin = 1, MAXGEN = 1000, MAXSIZE = 1000, NIND = 50, SUBPOP = 1, GGAP = 1, selectStyle = 'tour', recombinStyle = 'xovdprs', recopt = 0.9, pm = None, drawing = 1)

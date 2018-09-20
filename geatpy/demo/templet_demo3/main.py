import numpy as np
import geatpy as ga

# 获取函数接口地址
AIM_M = __import__('aimfuc')
# 变量设置
ranges = np.vstack([np.zeros((1, 4)), np.ones((1, 4))]) # 生成自变量的范围矩阵
borders = np.vstack([np.ones((1, 4)), np.ones((1, 4))]) # 生成自变量的边界矩阵
FieldDR = ga.crtfld(ranges, borders) # 生成区域描述器
# 调用编程模板
[pop_trace, var_trace, times] = ga.sga_real_templet(AIM_M, 'aimfuc', None, None, FieldDR, problem = 'I', maxormin = -1, MAXGEN = 50, NIND = 10, SUBPOP = 1, GGAP = 0.9, selectStyle = 'sus', recombinStyle = 'xovdp', recopt = 0.9, pm = 0.1)

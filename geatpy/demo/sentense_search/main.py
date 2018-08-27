# -*- coding: utf-8 -*-
"""             *****main.py*****
use GA algorithm to search the sentence: Tom is a little boy.
"""
import numpy as np
import geatpy as ga

strs = 'Tom is a little boy.' # define the sentence
words = []
for c in strs:
    words.append(ord(c)) # change the words to ascii code.

def aim(Phen): # define the aim function
    real = words
    diff = np.sum((Phen - real)**2, 1)
    return np.array([diff]).T

if __name__ == "__main__":
    AIM_M = __import__('main') # 获取函数接口所在文件的地址
    # 变量设置
    ranges = np.vstack([32 * np.ones((1, len(words))), 122 * np.ones((1, len(words)))]) # 生成自变量的范围矩阵
    borders = np.vstack([np.ones((1, len(words))), 122 * np.ones((1, len(words)))]) # 生成自变量的边界矩阵
    FieldDR = ga.crtfld(ranges, borders) # 生成区域描述器
    # 调用编程模板
    [pop_trace, var_trace, times] = ga.sga_new_real_templet(AIM_M, 'aim', None, None, FieldDR, problem = 'I', maxormin = 1, MAXGEN = 2000, NIND = 50, SUBPOP = 1, GGAP = 0.9, selectStyle = 'sus', recombinStyle = 'xovdp', recopt = 0.9, pm = 0.1, drawing = 1)
    # 输出结果
    for num in var_trace[np.argmin(pop_trace[:, 1]), :]:
        print(chr(int(num)), end = '')
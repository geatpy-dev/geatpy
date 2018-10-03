# -*- coding: utf-8 -*-
"""
执行脚本quickstart_demo.py
描述:
    本demo展示了如何不使用Geatpy提供的进化算法框架、纯粹地写编程脚本来解决一个无约束的单目标优化问题。
"""

import numpy as np
import geatpy as ga # 导入geatpy库
import time

"""============================函数定义============================"""
def aimfuc(Phen): # 传入种群染色体矩阵解码后的基因表现型矩阵
    x = np.array([Phen[:, 0]]).T # 取出Phen第一列得到种群所有个体的x值
    y = np.array([Phen[:, 1]]).T # 取出Phen第二列得到种群所有个体的y值
    # 计算MoCormick函数值
    f = np.sin(x + y) + (x - y) ** 2 - 1.5 * x + 2.5 * y +1
    return f
"""============================变量设置============================"""
x1 = [-1.5, 4]; x2 = [-3, 4]       # 自变量范围
b1 = [1, 1];    b2 = [1, 1]        # 自变量边界
codes = [1, 1]                     # 变量的编码方式，2个变量均使用格雷编码
precisions =[6, 6]                 # 变量的精度，10表示精确到小数点后10位
scales = [0, 0]                    # 采用算术刻度
ranges=np.vstack([x1, x2]).T       # 生成自变量的范围矩阵
borders=np.vstack([b1, b2]).T      # 生成自变量的边界矩阵
"""========================遗传算法参数设置========================="""
NIND = 50                # 种群规模
MAXGEN = 1000            # 最大遗传代数
GGAP = 0.8               # 代沟：子代与父代个体不相同的概率为0.8
selectStyle = 'sus';     # 遗传算法的选择方式设为"sus"——随机抽样选择
recombinStyle = 'xovdp'  # 遗传算法的重组方式，设为两点交叉
recopt = 0.9             # 交叉概率                                                      
pm = 0.1                 # 变异概率
SUBPOP = 1               # 设置种群数为1
"""=========================开始遗传算法进化========================"""
FieldD = ga.crtfld(ranges,borders,precisions,codes,scales) # 调用函数创建区域描述器
Lind = np.sum(FieldD[0, :]) # 计算编码后的染色体长度
Chrom = ga.crtbp(NIND, Lind) # 根据区域描述器生成二进制种群
Phen = ga.bs2rv(Chrom, FieldD) #对初始种群进行解码
ObjV = aimfuc(Phen) # 计算初始种群个体的目标函数值
# 定义进化记录器，初始值为nan
pop_trace = (np.zeros((MAXGEN, 2)) * np.nan)
# 定义种群最优个体记录器，记录每一代最优个体的染色体，初始值为nan
ind_trace = (np.zeros((MAXGEN, Lind)) * np.nan)
# 开始进化！！
start_time = time.time() # 开始计时
for gen in range(MAXGEN):
    FitnV = ga.ranking(ObjV) # 根据目标函数大小分配适应度值
    SelCh=ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP) # 选择
    SelCh=ga.recombin(recombinStyle, SelCh, recopt, SUBPOP) #交叉
    SelCh=ga.mutbin(SelCh, pm) # 二进制种群变异
    Phen = ga.bs2rv(SelCh, FieldD) # 对育种种群进行解码(二进制转十进制)
    ObjVSel = aimfuc(Phen) # 求育种个体的目标函数值
    [Chrom,ObjV] = ga.reins(Chrom,SelCh,SUBPOP,2,1,ObjV,ObjVSel) # 重插入得到新一代种群
    # 记录
    best_ind = np.argmin(ObjV) # 计算当代最优个体的序号
    pop_trace[gen, 0] = ObjV[best_ind] # 记录当代种群最优个体目标函数值
    pop_trace[gen, 1] = np.sum(ObjV) / ObjV.shape[0] # 记录当代种群的目标函数均值
    ind_trace[gen, :] = Chrom[best_ind, :] # 记录当代种群最优个体的变量值
# 进化完成
end_time = time.time() # 结束计时

"""============================绘图================================"""
ga.trcplot(pop_trace, [['最优个体目标函数值','种群的目标函数均值']], ['demo_result'])

"""============================输出结果============================"""
best_gen = np.argmin(pop_trace[:, 0]) # 计算最优种群是在哪一代
print('最优的目标函数值为：', np.min(pop_trace[:, 0]))
print('最优的控制变量值为：')
# 最优个体记录器存储的是各代种群最优个体的染色体，此处需要解码得到对应的基因表现型
variables = ga.bs2rv(ind_trace, FieldD) # 解码
for i in range(variables.shape[1]):
    print(variables[best_gen, i])
print('最优的一代是第',best_gen + 1,'代')
print('用时：', end_time - start_time, '秒')

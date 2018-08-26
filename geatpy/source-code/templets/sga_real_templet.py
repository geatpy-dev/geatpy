# -*- coding: utf-8 -*-

import numpy as np
import geatpy as ga
import time

def sga_real_templet(AIM_M, AIM_F, PUN_M, PUN_F, FieldDR, problem, maxormin, MAXGEN, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, drawing = 1):
    
    """
sga_real_templet.py - 单目标编程模板(实值编码)

语法：
    该函数除了参数drawing外，不设置可缺省参数。当某个参数需要缺省时，在调用函数时传入None即可。
    比如当没有罚函数时，则在调用编程模板时将第3、4个参数设置为None即可，如：
    sga_real_templet(AIM_M, 'aimfuc', None, None, ..., maxormin)

输入参数：
    AIM_M - 目标函数的地址，传入该函数前通常由AIM_M = __import__('目标函数名')语句得到
    
    AIM_F : str - 目标函数名
    
    PUN_M - 罚函数的地址，传入该函数前通常由PUN_M = __import__('罚函数名')语句得到
    
    PUN_F : str - 罚函数名
    
    FieldDR : array - 实际值种群区域描述器
        [lb;		(float) 指明每个变量使用的下界
         ub]		(float) 指明每个变量使用的上界
         注：不需要考虑是否包含变量的边界值。在crtfld中已经将是否包含边界值进行了处理
         本函数生成的矩阵的元素值在FieldDR的[下界, 上界)之间
    
    problem : str - 表明是整数问题还是实数问题，'I'表示是整数问题，'R'表示是实数问题                 
    
    maxormin int - 最小最大化标记，1表示目标函数最小化；-1表示目标函数最大化
    
    MAXGEN : int - 最大遗传代数
    
    NIND : int - 种群规模，即种群中包含多少个个体
    
    SUBPOP : int - 子种群数量，即对一个种群划分多少个子种群
    
    GGAP : float - 代沟，表示子代与父代染色体及性状不相同的概率
    
    selectStyle : str - 指代所采用的低级选择算子的名称，如'rws'(轮盘赌选择算子)
    
    recombinStyle: str - 指代所采用的低级重组算子的名称，如'xovsp'(单点交叉)
    
    recopt : float - 交叉概率
    
    pm : float - 重组概率
    
    drawing : int - (可选参数)，0表示不绘图，1表示绘制最终结果图。默认drawing为1

输出参数：
    pop_trace : array - 种群进化记录器(进化追踪器),
                        第0列记录着各代种群最优个体的目标函数值
                        第1列记录着各代种群的适应度均值
                        第2列记录着各代种群最优个体的适应度值
    
    var_trace : array - 变量记录器，记录着各代种群最优个体的变量值，每一列对应一个控制变量
    
    times     : float - 进化所用时间

"""
    """==========================初始化配置==========================="""
    # 获取目标函数和罚函数
    aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
    if PUN_F is not None:
        punishing = getattr(PUN_M, PUN_F) # 获得罚函数
    NVAR = FieldDR.shape[1] # 得到控制变量的个数
    # 定义进化记录器，初始值为nan
    pop_trace = (np.zeros((MAXGEN ,2)) * np.nan)
    # 定义变量记录器，记录控制变量值，初始值为nan
    var_trace = (np.zeros((MAXGEN ,NVAR)) * np.nan) 
    """=========================开始遗传算法进化======================="""
    if problem == 'R':
        Chrom = ga.crtrp(NIND, FieldDR) # 生成初始种群
    elif problem == 'I':
        Chrom = ga.crtip(NIND, FieldDR)
    ObjV = aimfuc(Chrom) # 求种群的目标函数值
    
    start_time = time.time() # 开始计时
    # 开始进化！！
    for gen in range(MAXGEN):
        FitnV = ga.ranking(maxormin * ObjV, None, SUBPOP)
        if PUN_F is not None:
            FitnV = punishing(Chrom, FitnV) # 调用罚函数
        # 进行遗传算子
        SelCh=ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP) # 选择
        SelCh=ga.recombin(recombinStyle, Chrom, recopt, SUBPOP) # 重组
        if problem == 'R':
            SelCh=ga.mutbga(SelCh,FieldDR, pm) # 变异
        elif problem == 'I':
            SelCh=ga.mutint(SelCh, FieldDR, pm)
        ObjVSel = aimfuc(SelCh) # 求后代的目标函数值
        # 重插入
        [Chrom, ObjV]=ga.reins(Chrom, SelCh, SUBPOP, 2, 1, ObjV, ObjVSel)
        pop_trace[gen,0] = np.sum(ObjV) / ObjV.shape[0] # 记录种群个体平均目标函数值
        if maxormin == 1:
            pop_trace[gen,1] = np.min(ObjV) # 记录当代目标函数的最优值
            var_trace[gen,:] = Chrom[np.argmin(ObjV), :] # 记录当代最优的控制变量值
        elif maxormin == -1:
            pop_trace[gen,1] = np.max(ObjV)
            var_trace[gen,:] = Chrom[np.argmax(ObjV), :] # 记录当代最优的控制变量值
    end_time = time.time() # 结束计时
    
    # 绘图
    if drawing == 1:
        ga.trcplot(pop_trace, [['种群个体平均目标函数值', '种群最优个体目标函数值']])
    # 输出结果
    if maxormin == 1:
        best_gen = np.argmin(pop_trace[:, 1]) # 记录最优种群是在哪一代
        print('最优的目标函数值为：', np.min(pop_trace[:, 1]))
    elif maxormin == -1:
        best_gen = np.argmax(pop_trace[:, 1]) # 记录最优种群是在哪一代
        print('最优的目标函数值为：', np.max(pop_trace[:, 1]))
    print('最优的控制变量值为：')
    for i in range(NVAR):
        print(var_trace[best_gen, i])
    print('最优的一代是第', best_gen + 1, '代')
    times = end_time - start_time
    print('时间已过', times, '秒')
    # 返回进化记录器、变量记录器以及执行时间
    return [pop_trace, var_trace, times]



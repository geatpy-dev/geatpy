# -*- coding: utf-8 -*-
"""自定义进化算法模板mintemp1.py"""
import numpy as np
import geatpy as ga # 导入geatpy库
import time

def mintemp1(AIM_M, AIM_F, PUN_M, PUN_F, ranges, borders, MAXGEN, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, maxormin):
    """==========================初始化配置==========================="""
    # 获取目标函数和罚函数
    aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
    punishing = getattr(PUN_M, PUN_F) # 获得罚函数
    FieldDR = ga.crtfld(ranges, borders) # 初始化区域描述器
    NVAR = ranges.shape[1] # 得到控制变量的个数
    # 定义进化记录器，初始值为nan
    pop_trace = (np.zeros((MAXGEN ,3)) * np.nan).astype('int64')
    # 定义变量记录器，记录控制变量值，初始值为nan
    var_trace = (np.zeros((MAXGEN ,NVAR)) * np.nan).astype('int64')    
    """=========================开始遗传算法进化======================="""
    Chrom = ga.crtip(NIND, FieldDR) # 根据区域描述器FieldDR生成整数型初始种群
    LegV = np.ones((NIND, 1)) # 生成可行性列向量，元素为1表示对应个体是可行解，0表示非可行解
    [ObjV, LegV] = aimfuc(Chrom, LegV) # 计算种群目标函数值，同时更新LegV
    start_time = time.time() # 开始计时
    # 开始进化！！
    for gen in range(MAXGEN):
        FitnV = ga.ranking(maxormin * ObjV, LegV) # 计算种群适应度
        FitnV = punishing(LegV, FitnV) # 调用罚函数
        # 记录进化过程
        bestIdx = np.argmax(FitnV)
        if LegV[bestIdx] != 0:
            feasible = np.where(LegV != 0)[0] # 排除非可行解
            # 记录当代种群的适应度均值
            pop_trace[gen, 1] = np.sum(FitnV[feasible]) / FitnV[feasible].shape[0]
            # 记录当代种群最优个体的目标函数值
            pop_trace[gen, 0] = ObjV[bestIdx]
            # 记录当代种群的最优个体的适应度值
            pop_trace[gen, 2] = FitnV[bestIdx]
            # 记录当代种群最优个体的变量值
            var_trace[gen, :] = Chrom[bestIdx, :] 
        # 进行遗传操作！！
        SelCh=ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP) # 选择
        SelCh=ga.recombin(recombinStyle, SelCh, recopt, SUBPOP) #交叉
        SelCh=ga.mutint(SelCh, FieldDR, pm) # 实值变异
        LegVSel = np.ones((SelCh.shape[0], 1)) # 创建育种个体的可行性列向量
        [ObjVSel, LegVSel] = aimfuc(SelCh, LegVSel) # 求育种个体的目标函数值
        FitnVSel = punishing(LegVSel, FitnV) # 调用罚函数
        [Chrom,ObjV,LegV] = ga.reins(Chrom,SelCh,SUBPOP,1,1,FitnV,FitnVSel,ObjV,ObjVSel,
        LegV,LegVSel) #重插入
    end_time = time.time() # 结束计时
    # 返回进化记录器、变量记录器以及执行时间
    return [pop_trace, var_trace, end_time - start_time]
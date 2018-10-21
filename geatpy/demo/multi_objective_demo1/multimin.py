# -*- coding: utf-8 -*-
"""自定义的进化算法模板multimin.py"""
import numpy as np
import geatpy as ga # 导入geatpy库
import time

def multimin(AIM_M, AIM_F, PUN_M, PUN_F, NIND, NVAR, Base, MAXGEN, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, maxormin):
    """==========================初始化配置==========================="""
    # 获取目标函数和罚函数
    aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
    if PUN_F is not None:
        punishing = getattr(PUN_M, PUN_F) # 获得罚函数
    BaseV = ga.crtbase(NVAR, Base)
    """=========================开始遗传算法进化======================="""
    Chrom = ga.crtbp(NIND, BaseV) # 创建简单离散种群
    LegV = np.ones((NIND, 1)) # 初始化种群个体的可行性列向量
    [ObjV,LegV] = aimfuc(Chrom,LegV) # 计算种群目标函数值
    NDSet = np.zeros((0, Chrom.shape[1])) # 定义帕累托最优解集合(初始为空集)
    NDSetObjV = np.zeros((0, ObjV.shape[1])) # 定义帕累托最优解的目标函数值记录器
    start_time = time.time() # 开始计时
    # 开始进化！！
    for gen in range(MAXGEN):
        # 求种群的非支配个体以及基于被支配数的适应度
        [FitnV, frontIdx] = ga.ndominfast(maxormin * ObjV, LegV)
        if PUN_F is not None:
            FitnV = punishing(LegV, FitnV) # 调用罚函数，不满足约束条件的个体适应度被设为0
        # 更新帕累托最优集以及种群非支配个体的适应度
        [FitnV, NDSet, NDSetObjV, repnum] = ga.upNDSet(Chrom, maxormin * ObjV, FitnV, NDSet, maxormin * NDSetObjV, frontIdx, LegV)
        NDSetObjV *= maxormin # 还原在传入upNDSet函数前被最小化处理过的NDSetObjV
        # 进行遗传操作！！
        SelCh=ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP) # 选择
        SelCh=ga.recombin(recombinStyle, SelCh, recopt, SUBPOP) #交叉
        SelCh=ga.mut(SelCh, BaseV, pm) # 变异
        LegVSel = np.ones((SelCh.shape[0], 1)) # 初始化育种种群的可行性列向量
        [ObjVSel,LegVSel] = aimfuc(SelCh, LegVSel) # 求育种个体的目标函数值
        # 求种群的非支配个体以及基于被支配数的适应度
        [FitnVSel, frontIdx] = ga.ndominfast(maxormin * ObjVSel, LegVSel)
        if PUN_F is not None:
            FitnVSel = punishing(LegVSel, FitnVSel) # 调用罚函数
        [Chrom,ObjV,LegV] = ga.reins(Chrom,SelCh,SUBPOP,1,1,FitnV,FitnVSel,ObjV,ObjVSel,LegV,LegVSel) #重插入
    end_time = time.time() # 结束计时    
    # 返回进化记录器、变量记录器以及执行时间
    return [ObjV, NDSet, NDSetObjV, end_time - start_time]
"""自定义进化算法模板multimin.py"""
import numpy as np
import geatpy as ga # 导入geatpy库
import time

def multimin(AIM_M, AIM_F, NIND, NVAR, Base, MAXGEN, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, maxormin):
    """==========================初始化配置==========================="""
    # 获取目标函数和罚函数
    aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
    BaseV = ga.crtbase(NVAR, Base)
    """=========================开始遗传算法进化======================="""
    Chrom = ga.crtbp(NIND, BaseV) # 创建简单离散种群
    ObjV = aimfuc(Chrom) # 计算种群目标函数值
    NDSet = np.zeros((0, Chrom.shape[1])) # 定义帕累托最优解集合(初始为空集)
    NDSetObjV = np.zeros((0, ObjV.shape[1])) # 定义帕累托最优解的目标函数值记录器
    start_time = time.time() # 开始计时
    # 开始进化！！
    for gen in range(MAXGEN):
        # 求种群的非支配个体以及基于被支配数的适应度
        [FitnV, frontIdx] = ga.ndominfast(maxormin * ObjV)
        # 更新帕累托最优集以及种群非支配个体的适应度
        [FitnV, NDSet, NDSetObjV, repnum] = ga.upNDSet(Chrom, maxormin * ObjV, FitnV, NDSet, maxormin * NDSetObjV, frontIdx)
        # 进行遗传操作！！
        SelCh=ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP) # 选择
        SelCh=ga.recombin(recombinStyle, SelCh, recopt, SUBPOP) #交叉
        SelCh=ga.mut(SelCh, BaseV, pm) # 变异
        ObjVSel = aimfuc(SelCh) # 求育种个体的目标函数值
        # 求种群的非支配个体以及基于被支配数的适应度
        [FitnVSel, frontIdx] = ga.ndominfast(maxormin * ObjVSel)
        [Chrom,ObjV] = ga.reins(Chrom,SelCh,SUBPOP,1,1,FitnV,FitnVSel,ObjV,ObjVSel) #重插入
    end_time = time.time() # 结束计时    
    # 返回进化记录器、变量记录器以及执行时间
    return [ObjV, NDSet, NDSetObjV, end_time - start_time]
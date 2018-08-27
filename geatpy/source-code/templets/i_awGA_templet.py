# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ga # 导入geatpy库
import time

def i_awGA_templet(AIM_M, AIM_F, PUN_M, PUN_F, FieldDR, problem, maxormin, MAXGEN, MAXSIZE, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, drawing = 1):
    
    """
i_awGA_templet.py - 基于交互式适应性权重聚合法(i-awGA)求解多目标优化问题的编程模板

语法：
    该函数除参数drawing外，不设置可缺省参数。当某个参数需要缺省时，在调用函数时传入None即可。
    比如当没有罚函数时，则在调用编程模板时将第3、4个参数设置为None即可，如：
    i-awGA_templet(AIM_M, 'aimfuc', None, None, ..., maxormin)
    
    本模板实现了基于改进的交互式适应性权重聚合法(i-awGA)的多目标优化搜索，
    改进之处是将原i-awGA中使用的Deb非支配排序法变为快速非支配排序法，
    并维护一个全局帕累托最优集来实现帕累托前沿的搜索，故并不需要保证种群所有个体都是非支配的
    
输入参数：
    AIM_M - 目标函数的地址，传入该函数前通常由AIM_M = __import__('目标函数名')语句得到
            目标函数规范定义：f = aimfuc(Phen)
            其中Phen是种群的表现型矩阵
    
    AIM_F : str - 目标函数名
    
    PUN_M - 罚函数的地址，传入该函数前通常由PUN_M = __import__('罚函数名')语句得到
            罚函数规范定义： f = punishing(Phen, FitnV)
            其中Phen是种群的表现型矩阵, FitnV为种群个体适应度列向量
    
    PUN_F : str - 罚函数名
    
    FieldDR : array - 实际值种群区域描述器
        [lb;		(float) 指明每个变量使用的下界
         ub]		(float) 指明每个变量使用的上界
         注：不需要考虑是否包含变量的边界值。在crtfld中已经将是否包含边界值进行了处理
         本函数生成的矩阵的元素值在FieldDR的[下界, 上界)之间
    
    problem : str - 表明是整数问题还是实数问题，'I'表示是整数问题，'R'表示是实数问题
    
    maxormin int - 最小最大化标记，1表示目标函数最小化；-1表示目标函数最大化
    
    MAXGEN : int - 最大遗传代数
    
    MAXSIZE : int - 帕累托最优集最大规模
    
    NIND : int - 种群规模，即种群中包含多少个个体
    
    SUBPOP : int - 子种群数量，即对一个种群划分多少个子种群
    
    GGAP : float - 代沟，表示子代与父代染色体及性状不相同的概率
    
    selectStyle : str - 指代所采用的低级选择算子的名称，如'rws'(轮盘赌选择算子)
    
    recombinStyle: str - 指代所采用的低级重组算子的名称，如'xovsp'(单点交叉)
    
    recopt : float - 交叉概率
    
    pm : float - 重组概率
    
    drawing : int - (可选参数)，0表示不绘图，1表示绘制最终结果图，2表示绘制进化过程的动画。
                    默认drawing为1

算法描述:
    本模板实现了基于适应性权重聚合法(awGA)的多目标优化搜索，
    通过维护一个全局帕累托最优集来实现帕累托前沿的搜索，故并不需要保证种群所有个体都是非支配的
    
"""
    
    #==========================初始化配置===========================
    # 获取目标函数和罚函数
    aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
    if PUN_F is not None:
        punishing = getattr(PUN_M, PUN_F) # 获得罚函数
    #=========================开始遗传算法进化=======================
    if problem == 'R':
        Chrom = ga.crtrp(NIND, FieldDR) # 生成实数值种群
    elif problem == 'I':
        Chrom = ga.crtip(NIND, FieldDR) # 生成整数值种群
    ObjV = aimfuc(Chrom) # 计算种群目标函数值
    NDSet = np.zeros((0, Chrom.shape[1])) # 定义帕累托最优解记录器
    NDSetObjV = np.zeros((0, ObjV.shape[1])) # 定义帕累托最优解的目标函数值记录器
    ax = None
    start_time = time.time() # 开始计时
    # 开始进化！！
    for gen in range(MAXGEN):
        if NDSet.shape[0] > MAXSIZE:
            break
        [CombinObjV, weight] = ga.awGA(ObjV) # 适应性权重法求聚合目标函数值
        FitnV  = ga.ranking(maxormin * CombinObjV) # 根据加权单目标计算适应度
        if PUN_F is not None:
            FitnV = punishing(Chrom, FitnV) # 调用罚函数
        [FitnV, frontIdx] = ga.ndominfast(maxormin * ObjV, FitnV) # 求种群的非支配个体，并更新适应度
        # 更新帕累托最优集以及种群非支配个体的适应度
        [FitnV, NDSet, NDSetObjV, repnum] = ga.upNDSet(Chrom, maxormin * ObjV, FitnV, NDSet, maxormin * NDSetObjV, frontIdx)
        # 进行遗传操作！！
        SelCh=ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP) # 选择
        SelCh=ga.recombin(recombinStyle, SelCh, recopt, SUBPOP) #交叉
        if problem == 'R':
            SelCh=ga.mutbga(SelCh,FieldDR, pm) # 变异
            if repnum > Chrom.shape[0] * 0.05: # 当最优个体重复率高达5%时，进行一次高斯变异
                SelCh=ga.mutgau(SelCh, FieldDR, pm) # 高斯变异
        elif problem == 'I':
            SelCh=ga.mutint(SelCh, FieldDR, pm)
        ObjVSel = aimfuc(SelCh) # 求育种个体的目标函数值
        [CombinObjV, weight] = ga.awGA(maxormin * ObjVSel) # 适应性权重法求聚合目标函数值
        FitnVSel = ga.ranking(maxormin * CombinObjV) # 根据聚合目标求育种个体适应度
        if PUN_F is not None:
            FitnVSel = punishing(Chrom, FitnVSel) # 调用罚函数
        [Chrom,ObjV] = ga.reins(Chrom,SelCh,SUBPOP,1,0.9,FitnV,FitnVSel,ObjV,ObjVSel) #重插入
        if drawing == 2:
            ax = ga.frontplot(NDSetObjV, False, ax, gen + 1) # 绘制动态图
    end_time = time.time() # 结束计时
    #=========================绘图及输出结果=========================
    if drawing != 0:
        ga.frontplot(NDSetObjV,True)
    times = end_time - start_time
    print('用时：', times, '秒')
    print('帕累托前沿点个数：', NDSet.shape[0], '个')
    print('单位时间找到帕累托前沿点个数：', int(NDSet.shape[0] // times), '个')
    # 返回帕累托最优集以及执行时间
    return [ObjV, NDSet, NDSetObjV, end_time - start_time]
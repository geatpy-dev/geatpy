# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ga # 导入geatpy库
import time

def nsga2_templet(AIM_M, AIM_F, PUN_M, PUN_F, FieldDR, problem, maxormin, MAXGEN, MAXSIZE, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, drawing = 1):
    
    """
nsga2_templet.py - 基于改进NSGA-Ⅱ算法求解多目标优化问题编程模板

语法：
    该函数除参数drawing外，不设置可缺省参数。当某个参数需要缺省时，在调用函数时传入None即可。
    比如当没有罚函数时，则在调用编程模板时将第3、4个参数设置为None即可，如：
    nsga2_templet(AIM_M, 'aimfuc', None, None, ..., maxormin)
    
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
    
    GGAP : float - 代沟，本模板中该参数为无用参数，仅为了兼容同类模板而设
    
    selectStyle : str - 指代所采用的低级选择算子的名称，如'rws'(轮盘赌选择算子)
    
    recombinStyle: str - 指代所采用的低级重组算子的名称，如'xovsp'(单点交叉)
    
    recopt : float - 交叉概率
    
    pm : float - 重组概率
    
    drawing : int - (可选参数)，0表示不绘图，1表示绘制最终结果图，2表示绘制进化过程的动画。
                    默认drawing为1
算法描述:
    传统NSGA-Ⅱ算法的帕累托最优解来只源于当代种群个体，这样难以高效地获取更多的帕累托最优解
    同时难以把种群大小控制在合适的范围内，
    改进的NSGA2整体上沿用传统的NSGA-Ⅱ算法，
    不同的是，该算法通过维护一个全局帕累托最优集来实现帕累托前沿的搜索，
    故并不需要保证种群所有个体都是非支配的。
    
值得注意的是：
    尽管当全局帕累托最优集大小比种群规模大时，算法的时间复杂度比原NSGA-Ⅱ算法要高，
    但算法整体的时间复杂度要低与原NSGA-Ⅱ算法，
    因此单位时间内生成的无重复帕累托最优解个数要多于原NSGA-Ⅱ算法
    
"""
    # 获取目标函数和罚函数
    aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
    if PUN_F is not None:
        punishing = getattr(PUN_M, PUN_F) # 获得罚函数
    #==========================初始化配置===========================
    GGAP = 0.5 # 为了避免父子两代合并后种群数量爆炸，要让代沟为0.5
    # 获取目标函数和罚函数
    aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
    #=========================开始遗传算法进化=======================
    if problem == 'R':
        Chrom = ga.crtrp(NIND, FieldDR) # 生成实数值种群
    elif problem == 'I':
        Chrom = ga.crtip(NIND, FieldDR) # 生成整数值种群
    ObjV = aimfuc(Chrom) # 计算种群目标函数值
    NDSet = np.zeros((0, Chrom.shape[1])) # 定义帕累托最优解集合(初始为空集)
    NDSetObjV = np.zeros((0, ObjV.shape[1])) # 定义帕累托最优解对应的目标函数集合(初始为空集)
    ax = None
    start_time = time.time() # 开始计时
    [FitnV, levels] = ga.ndomindeb(maxormin * ObjV, 1) # deb非支配分级
    if PUN_F is not None:
        FitnV = punishing(Chrom, FitnV) # 调用罚函数
    frontIdx = np.where(levels == 1)[0] # 处在第一级的个体即为种群的非支配个体
        # 更新帕累托最优集以及种群非支配个体的适应度
    [FitnV, NDSet, NDSetObjV, repnum] = ga.upNDSet(Chrom, maxormin * ObjV, FitnV, NDSet, maxormin * NDSetObjV, frontIdx)
    # 开始进化！！
    for gen in range(MAXGEN):
        if NDSet.shape[0] > MAXSIZE:
            break
        # 进行遗传操作！！
        SelCh=ga.recombin(recombinStyle, Chrom, recopt, SUBPOP) #交叉
        if problem == 'R':
            SelCh=ga.mutbga(SelCh,FieldDR, pm) # 变异
            if repnum > Chrom.shape[0] * 0.05: # 当最优个体重复率高达5%时，进行一次高斯变异
                SelCh=ga.mutgau(SelCh, FieldDR, pm) # 高斯变异
        elif problem == 'I':
            SelCh=ga.mutint(SelCh, FieldDR, pm)
        # 父子合并
        Chrom = np.vstack([Chrom, SelCh])
        ObjV =  aimfuc(Chrom) # 求目标函数值
        [FitnV, levels] = ga.ndomindeb(maxormin * ObjV, 1) # deb非支配分级
        if PUN_F is not None:
            FitnV = punishing(Chrom, FitnV) # 调用罚函数
        frontIdx = np.where(levels == 1)[0] # 处在第一级的个体即为种群的非支配个体
        # 更新帕累托最优集以及种群非支配个体的适应度
        [FitnV, NDSet, NDSetObjV, repnum] = ga.upNDSet(Chrom, maxormin * ObjV, FitnV, NDSet, maxormin * NDSetObjV, frontIdx)
        # 计算每个目标下个体的聚集距离(不需要严格计算欧氏距离，计算绝对值即可)
        for i in range(ObjV.shape[1]):
            idx = np.argsort(ObjV[:, i], 0)
            dis = np.abs(np.diff(ObjV[idx, i].T, 1).T) / (np.max(ObjV[idx, i]) - np.min(ObjV[idx, i]) + 1) # 差分计算距离
            dis = np.hstack([dis, dis[-1]])
            FitnV[idx, 0] += dis # 根据聚集距离修改适应度，以增加种群的多样性
        Chrom=ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP) # 选择出下一代
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
    return [ObjV, NDSet, NDSetObjV, times]

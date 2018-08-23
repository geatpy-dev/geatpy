# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ga # 导入geatpy库
import time

def nsga2_templet(AIM_M, AIM_F, PUN_M, PUN_F, ranges, borders, precisions, maxormin, MAXGEN, MAXSIZE, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, drawing = 1):
    
    """
nsga2_templet.py - 基于改进NSGA-Ⅱ算法求解多目标优化问题编程模板

语法：
    该函数除参数drawing外，不设置可缺省参数。当某个参数需要缺省时，在调用函数时传入None即可。
    比如当没有罚函数时，则在调用编程模板时将第3、4个参数设置为None即可，如：
    nsga2_templet(AIM_M, 'aimfuc', None, None, ..., maxormin)
    
输入参数：
    AIM_M - 目标函数的地址，传入该函数前通常由AIM_M = __import__('目标函数名')语句得到
    
    AIM_F : str - 目标函数名
    
    PUN_M - 罚函数的地址，传入该函数前通常由PUN_M = __import__('罚函数名')语句得到
    
    PUN_F : str - 罚函数名
    
    ranges : array  - 代表自变量的范围矩阵，要求上界必须大于下界
        例如：[[1, 2, 3],
              [3, 4, 5]]
        表示有3个控制变量，其范围分别是1-3, 2-4, 3-5
                         
    borders : list -（可选参数）代表是否包含变量范围的边界，为1代表控制变量的范围包含该边界
        当为None时，默认设置为全是1的矩阵
        例如：[[1, 0, 1],
              [0, 1, 1]]
        表示上面的三个控制变量的范围分别是：[1, 3)、(2, 4]、[3, 5]
    
    precisions : list -（可选参数）代表控制变量的精度，
        如等于4，表示对应的控制变量的编码可以精确到小数点后4位。
        当precisions为None时，默认precision为1*n的0矩阵(此时表示种群是离散编码的)
        precision的元素必须不小于0
    
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
    FieldDR = ga.crtfld(ranges, borders, precisions)
    #=========================开始遗传算法进化=======================
    Chrom = ga.crtrp(NIND, FieldDR) # 创建简单离散种群
    ObjV = aimfuc(Chrom) # 计算种群目标函数值
    NDSet = np.zeros((0, ObjV.shape[1])) # 定义帕累托最优解集合(初始为空集)
    ax = None
    start_time = time.time() # 开始计时
    [FitnV, levels] = ga.ndomindeb(maxormin * ObjV, 1) # deb非支配分级
    if PUN_F is not None:
        FitnV = punishing(Chrom, FitnV) # 调用罚函数
    frontIdx = np.where(levels == 1)[0] # 处在第一级的个体即为种群的非支配个体
        # 更新帕累托最优集以及种群非支配个体的适应度
    [FitnV, NDSet, repnum] = ga.upNDSet(FitnV, maxormin * ObjV, maxormin * NDSet, frontIdx)
    # 开始进化！！
    for gen in range(MAXGEN):
        if NDSet.shape[0] > MAXSIZE:
            break
        # 进行遗传操作！！
        SelCh=ga.recombin(recombinStyle, Chrom, recopt, SUBPOP) #交叉
        SelCh=ga.mutbga(SelCh, FieldDR, pm) # 变异
        if repnum > Chrom.shape[0] * 0.05: # 当最优个体重复率高达5%时，进行一次高斯变异
            SelCh=ga.mutgau(SelCh, FieldDR, pm) # 高斯变异
        # 父子合并
        Chrom = np.vstack([Chrom, SelCh])
        ObjV =  aimfuc(Chrom) # 求目标函数值
        [FitnV, levels] = ga.ndomindeb(maxormin * ObjV, 1) # deb非支配分级
        if PUN_F is not None:
            FitnV = punishing(Chrom, FitnV) # 调用罚函数
        frontIdx = np.where(levels == 1)[0] # 处在第一级的个体即为种群的非支配个体
        # 更新帕累托最优集以及种群非支配个体的适应度
        [FitnV, NDSet, repnum] = ga.upNDSet(FitnV, maxormin * ObjV, maxormin * NDSet, frontIdx)
        # 计算每个目标下个体的聚集距离(不需要严格计算欧氏距离，计算绝对值即可)
        for i in range(ObjV.shape[1]):
            idx = np.argsort(ObjV[:, i], 0)
            dis = np.abs(np.diff(ObjV[idx, i].T, 1).T) / (np.max(ObjV[idx, i]) - np.min(ObjV[idx, i]) + 1) # 差分计算距离
            dis = np.hstack([dis, dis[-1]])
            FitnV[idx, 0] += dis # 根据聚集距离修改适应度，以增加种群的多样性
        Chrom=ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP) # 选择出下一代
        if drawing == 2:
            ax = ga.frontplot(NDSet, False, ax, gen + 1) # 绘制动态图
    end_time = time.time() # 结束计时
    #=========================绘图及输出结果=========================
    if drawing != 0:
        ga.frontplot(NDSet,True)
    times = end_time - start_time
    print('用时：' + str(times) + '秒')
    print('帕累托前沿点个数：' + str(NDSet.shape[0]) + '个')
    print('单位时间找到帕累托前沿点个数：' + str(NDSet.shape[0] // times) + '个')
    # 返回帕累托最优集以及执行时间
    return [ObjV, NDSet, times]

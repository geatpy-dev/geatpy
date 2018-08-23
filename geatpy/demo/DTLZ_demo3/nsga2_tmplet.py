"""
nsga2_templet.py - 基于改进NSGA-Ⅱ算法求解多目标优化问题编程模板

算法描述:
    传统NSGA-Ⅱ算法的帕累托最优解来只源于当代种群个体，这样难以高效地获取更多的帕累托最优解
    同时难以把种群大小控制在合适的范围内，
    改进的NSGA2整体上沿用传统的NSGA-Ⅱ算法，
    不同的是，该算法通过维护一个全局帕累托最优集来实现帕累托前沿的搜索，
    故并不需要保证种群所有个体都是非支配的。
    算法的时间复杂度为O(r*N*max(N,S)),其中S为全局帕累托最优集的大小,N为种群规模,r为目标数。

值得注意的是：
    尽管当全局帕累托最优集大小比种群规模大时，算法的时间复杂度比原NSGA-Ⅱ算法要高，
    但算法整体的时间复杂度要低与原NSGA-Ⅱ算法，
    因此单位时间内生成的无重复帕累托最优解个数要多于原NSGA-Ⅱ算法
    
"""
import numpy as np
import geatpy as ga # 导入geatpy库
import time

def nsga2(AIM_M, AIM_F, NIND, ranges, borders, precisions, MAXGEN, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, maxormin):
    """==========================初始化配置==========================="""
    GGAP = 0.5 * GGAP # 为了避免父子两代合并后种群数量爆炸，要让代沟不超过0.5
    # 获取目标函数和罚函数
    aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
    FieldDR = ga.crtfld(ranges, borders, precisions)
    """=========================开始遗传算法进化======================="""
    Chrom = ga.crtrp(NIND, FieldDR) # 创建简单离散种群
    ObjV = aimfuc(Chrom) # 计算种群目标函数值
    NDSet = np.zeros((0, ObjV.shape[1])) # 定义帕累托最优解集合(初始为空集)
    
    start_time = time.time() # 开始计时
    
    [FitnV, levels, maxLevel] = ga.ndomindeb(maxormin * ObjV, 1) # deb非支配分级
    frontIdx = np.where(levels == 1)[0] # 处在第一级的个体即为种群的非支配个体
        # 更新帕累托最优集以及种群非支配个体的适应度
    [FitnV, NDSet, repnum] = ga.upNDSet(FitnV, maxormin * ObjV, maxormin * NDSet, frontIdx)
    
    # 开始进化！！
    for gen in range(MAXGEN):
        if NDSet.shape[0] > 2 * ObjV.shape[0]:
            break
        # 进行遗传操作！！
        SelCh=ga.recombin(recombinStyle, Chrom, recopt, SUBPOP) #交叉
        SelCh=ga.mutbga(SelCh, FieldDR, pm) # 变异
        if repnum > Chrom.shape[0] * 0.05: # 当最优个体重复率高达5%时，进行一次高斯变异
            SelCh=ga.mutgau(SelCh, FieldDR, pm) # 高斯变异
        # 父子合并
        Chrom = np.vstack([Chrom, SelCh])
        ObjV =  aimfuc(Chrom) # 求目标函数值
        [FitnV, levels, maxLevel] = ga.ndomindeb(maxormin * ObjV, 1) # deb非支配分级
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
        
    end_time = time.time() # 结束计时   
    
    # 返回帕累托最优集以及执行时间
    return [ObjV, NDSet, end_time - start_time]

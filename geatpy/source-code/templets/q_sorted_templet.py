# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ga # 导入geatpy库
import time

def q_sorted_templet(AIM_M, AIM_F, PUN_M, PUN_F, ranges, borders, precisions, maxormin, MAXGEN, MAXSIZE, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, drawing = 1):
    
    """
q_sorted_templet.py - 基于交互式适应性权重聚合法(i-awGA)求解多目标优化问题的编程模板

语法：
    该函数除参数drawing外，不设置可缺省参数。当某个参数需要缺省时，在调用函数时传入None即可。
    比如当没有罚函数时，则在调用编程模板时将第3、4个参数设置为None即可，如：
    q_sorted_templet(AIM_M, 'aimfuc', None, None, ..., maxormin)
    
    本模板实现了基于改进的交互式适应性权重聚合法(i-awGA)的多目标优化搜索，
    改进之处是将原i-awGA中使用的Deb非支配排序法变为快速非支配排序法，
    并维护一个全局帕累托最优集来实现帕累托前沿的搜索，故并不需要保证种群所有个体都是非支配的
    
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
    
    GGAP : float - 代沟，表示子代与父代染色体及性状不相同的概率
    
    selectStyle : str - 指代所采用的低级选择算子的名称，如'rws'(轮盘赌选择算子)
    
    recombinStyle: str - 指代所采用的低级重组算子的名称，如'xovsp'(单点交叉)
    
    recopt : float - 交叉概率
    
    pm : float - 重组概率
    
    drawing : int - (可选参数)，0表示不绘图，1表示绘制最终结果图，2表示绘制进化过程的动画。
                    默认drawing为1
算法描述:
    本模板维护一个全局帕累托最优集来实现帕累托前沿的搜索
    利用快速非支配排序寻找每一代种群的非支配个体，并用它来不断更新全局帕累托最优集，
    故并不需要保证种群所有个体都是非支配的
    
"""
    
    #==========================初始化配置===========================
    # 获取目标函数和罚函数
    aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
    FieldDR = ga.crtfld(ranges, borders, precisions)
    #=========================开始遗传算法进化=======================
    Chrom = ga.crtrp(NIND, FieldDR) # 创建简单离散种群
    ObjV = aimfuc(Chrom) # 计算种群目标函数值
    NDSet = np.zeros((0, Chrom.shape[1])) # 定义帕累托最优解记录器
    NDSetObjV = np.zeros((0, ObjV.shape[1])) # 定义帕累托最优解的目标函数值记录器
    ax = None
    start_time = time.time() # 开始计时
    # 开始进化！！
    for gen in range(MAXGEN):
#        print(gen)
        if NDSet.shape[0] > MAXSIZE:
            break
        # 求种群的非支配个体以及基于被支配数的适应度
        [FitnV, frontIdx] = ga.ndominfast(maxormin * ObjV)
        # 更新帕累托最优集以及种群非支配个体的适应度
        [FitnV, NDSet, NDSetObjV, repnum] = ga.upNDSet(Chrom, maxormin * ObjV, FitnV, NDSet, maxormin * NDSetObjV, frontIdx)
#        if gen > 300:
#            return [Chrom, ObjV, FitnV, NDSet, NDSetObjV, frontIdx]
        # 进行遗传操作！！
        SelCh=ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP) # 选择
        SelCh=ga.recombin(recombinStyle, SelCh, recopt, SUBPOP) #交叉
        SelCh=ga.mutbga(SelCh, FieldDR, pm) # 变异
        if repnum > Chrom.shape[0] * 0.1: # 进行一次高斯变异
            SelCh=ga.mutgau(SelCh, FieldDR, pm) # 高斯变异
        ObjVSel = aimfuc(SelCh) # 求育种个体的目标函数值
        # 求种群的非支配个体以及基于被支配数的适应度
        [FitnVSel, frontIdx] = ga.ndominfast(maxormin * ObjVSel)
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
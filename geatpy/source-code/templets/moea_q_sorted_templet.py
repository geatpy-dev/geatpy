# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ga # 导入geatpy库
import time

def moea_q_sorted_templet(AIM_M, AIM_F, PUN_M, PUN_F, FieldDR, problem, maxormin, MAXGEN, MAXSIZE, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, distribute, drawing = 1):
    
    """
moea_q_sorted_templet.py - 基于快速非支配排序法求解多目标优化问题的进化算法模板

语法：
    该函数除参数drawing外，不设置可缺省参数。当某个参数需要缺省时，在调用函数时传入None即可。
    比如当没有罚函数时，则在调用编程模板时将第3、4个参数设置为None即可，如：
    moea_q_sorted_templet(AIM_M, 'aimfuc', None, None, ..., maxormin,...)

输入参数：
    AIM_M - 目标函数的地址，由AIM_M = __import__('目标函数所在文件名')语句得到
            目标函数规范定义：[f,LegV] = aimfuc(Phen,LegV)
            其中Phen是种群的表现型矩阵, LegV为种群的可行性列向量,f为种群的目标函数值矩阵
    
    AIM_F : str - 目标函数名
    
    PUN_M - 罚函数的地址，由PUN_M = __import__('罚函数所在文件名')语句得到
            罚函数规范定义： newFitnV = punishing(LegV, FitnV)
            其中LegV为种群的可行性列向量, FitnV为种群个体适应度列向量
            一般在罚函数中对LegV为0的个体进行适应度惩罚，返回修改后的适应度列向量newFitnV
    
    PUN_F : str - 罚函数名
    
    FieldDR : array - 实际值种群区域描述器
        [lb;		(float) 指明每个变量使用的下界
         ub]		(float) 指明每个变量使用的上界
         注：不需要考虑是否包含变量的边界值。在crtfld中已经将是否包含边界值进行了处理
         本函数生成的矩阵的元素值在FieldDR的[下界, 上界)之间
    
    problem : str - 表明是整数问题还是实数问题，'I'表示是整数问题，'R'表示是实数问题
    
    maxormin int - 最小最大化标记，1表示目标函数最小化；-1表示目标函数最大化
    
    MAXGEN : int - 最大遗传代数
    
    MAXSIZE : int - 帕累托最优集最大规模，当设为np.inf（无穷）时，模板不对帕累托最优解集规模作限制
    
    NIND : int - 种群规模，即种群中包含多少个个体
    
    SUBPOP : int - 子种群数量，即对一个种群划分多少个子种群
    
    GGAP : float - 代沟，表示子代与父代染色体及性状不相同的概率
    
    selectStyle : str - 指代所采用的低级选择算子的名称，如'rws'(轮盘赌选择算子)
    
    recombinStyle: str - 指代所采用的低级重组算子的名称，如'xovsp'(单点交叉)
    
    recopt : float - 交叉概率
    
    pm : float - 重组概率
    
    distribute : bool - 是否增强帕累托前沿的分布性（可能会造成收敛慢或帕累托前沿数目减少）
    
    drawing : int - (可选参数)，0表示不绘图，1表示绘制最终结果图，2表示绘制进化过程的动画。
                    默认drawing为1
算法描述:
    本模板维护一个全局帕累托最优集来实现帕累托前沿的搜索
    利用快速非支配排序寻找每一代种群的非支配个体，并用它来不断更新全局帕累托最优集，
    故并不需要保证种群所有个体都是非支配的

模板使用注意：
    1.本模板调用的目标函数形如：[ObjV,LegV] = aimfuc(Phen,LegV), 
      其中Phen表示种群的表现型矩阵, LegV为种群的可行性列向量(详见Geatpy数据结构)
    2.本模板调用的罚函数形如: newFitnV = punishing(LegV, FitnV), 
      其中FitnV为用其他算法求得的适应度
    若不符合上述规范，则请修改算法模板或自定义新算法模板
    3.关于'maxormin': geatpy的内核函数全是遵循“最小化目标”的约定的，即目标函数值越小越好。
      当需要优化最大化的目标时，需要设置'maxormin'为-1。
      本算法模板是正确使用'maxormin'的典型范例，其具体用法如下：
      当调用的函数传入参数包含与“目标函数值矩阵”有关的参数(如ObjV,ObjVSel,NDSetObjV等)时，
      查看该函数的参考资料(可用'help'命令查看，也可到官网上查看相应的教程)，
      里面若要求传入前对参数乘上'maxormin',则需要乘上。
      里面若要求对返回参数乘上'maxormin'进行还原，
      则调用函数返回得到的相应参数需要乘上'maxormin'进行还原，否则其正负号就会被改变。
    
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
    LegV = np.ones((NIND, 1)) # 初始化种群的可行性列向量
    [ObjV, LegV] = aimfuc(Chrom, LegV) # 计算种群目标函数值
    NDSet = np.zeros((0, Chrom.shape[1])) # 定义帕累托最优解记录器
    NDSetObjV = np.zeros((0, ObjV.shape[1])) # 定义帕累托最优解的目标函数值记录器
    ax = None # 存储上一帧动画
    start_time = time.time() # 开始计时
    # 开始进化！！
    for gen in range(MAXGEN):
        # 求种群的非支配个体以及基于被支配数的适应度
        [FitnV, frontIdx] = ga.ndominfast(maxormin * ObjV, LegV)
        if PUN_F is not None:
            FitnV = punishing(LegV, FitnV) # 调用罚函数作进一步的惩罚
        # 更新帕累托最优集以及种群非支配个体的适应度
        [FitnV, NDSet, NDSetObjV, repnum] = ga.upNDSet(Chrom, maxormin * ObjV, FitnV, NDSet, maxormin * NDSetObjV, frontIdx, LegV)
        NDSetObjV *= maxormin # 还原在传入upNDSet函数前被最小化处理过的NDSetObjV
        [NDSet, NDSetObjV] = ga.redisNDSet(NDSet, NDSetObjV, NDSetObjV.shape[1] * MAXSIZE) # 利用拥挤距离选择帕累托前沿的子集，在进化过程中最好比上限多筛选出几倍的点集
        if distribute == True: # 若要增强种群的分布性(可能会导致帕累托前沿搜索效率降低)
            # 计算每个目标下相邻个体的距离(不需要严格计算欧氏距离)
            for i in range(ObjV.shape[1]):
                idx = np.argsort(ObjV[:, i], 0)
                dis = np.diff(ObjV[idx, i]) / (np.max(ObjV[idx, i]) - np.min(ObjV[idx, i]) + 1) # 差分计算距离的偏移量占比，即偏移量除以目标函数的极差。加1是为了避免极差为0
                dis = np.hstack([dis, dis[-1]])
                FitnV[idx, 0] *= np.exp(dis) # 根据相邻距离修改适应度，突出相邻距离大的个体，以增加种群的多样性
        # 进行遗传操作！！
        SelCh=ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP) # 选择
        SelCh=ga.recombin(recombinStyle, SelCh, recopt, SUBPOP) #交叉
        if problem == 'R':
            SelCh=ga.mutbga(SelCh,FieldDR, pm) # 变异
            if repnum > Chrom.shape[0] * 0.01: # 当最优个体重复率高达1%时，进行一次高斯变异
                SelCh=ga.mutgau(SelCh, FieldDR, pm) # 高斯变异
        elif problem == 'I':
            SelCh=ga.mutint(SelCh, FieldDR, pm)
        LegVSel = np.ones((SelCh.shape[0], 1)) # 初始化育种种群的可行性列向量
        [ObjVSel, LegVSel] = aimfuc(SelCh, LegVSel) # 求育种个体的目标函数值
        # 求种群的非支配个体以及基于被支配数的适应度
        [FitnVSel, frontIdx] = ga.ndominfast(maxormin * ObjVSel, LegVSel)
        if PUN_F is not None:
            FitnVSel = punishing(LegVSel, FitnVSel) # 调用罚函数作进一步的惩罚
        [Chrom,ObjV,LegV] = ga.reins(Chrom,SelCh,SUBPOP,1,0.9,FitnV,FitnVSel,ObjV,ObjVSel,LegV,LegVSel) #重插入
        if drawing == 2:
            ax = ga.frontplot(NDSetObjV, False, ax, gen + 1) # 绘制动态图
    end_time = time.time() # 结束计时
    [NDSet, NDSetObjV] = ga.redisNDSet(NDSet, NDSetObjV, MAXSIZE) # 最后根据拥挤距离选择均匀分布的点
    #=========================绘图及输出结果=========================
    if drawing != 0:
        ga.frontplot(NDSetObjV,True)
    times = end_time - start_time
    print('用时：', times, '秒')
    print('帕累托前沿点个数：', NDSet.shape[0], '个')
    print('单位时间找到帕累托前沿点个数：', int(NDSet.shape[0] // times), '个')
    # 返回帕累托最优集以及执行时间
    return [ObjV, NDSet, NDSetObjV, end_time - start_time]

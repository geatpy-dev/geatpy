# -*- coding: utf-8 -*-

import numpy as np
import geatpy as ga
import time

def getUniBest(ChromBest):
    """
    该函数用于处理sol_trace，删除其中重复的解
    sol_trace是存储方程解的记录器，每一行对应一组解，每一列对应一组解里面的各个变量的值
    """
    
    if ChromBest.shape[0] > 1:
        ChromBest = np.sort(ChromBest, 0) # 对最优解进行排序
        dis = np.hstack([np.array([1]), np.sum(np.abs(np.diff(ChromBest.T, 1)), 0)]) # 计算相邻最优解之间的距离和
        dis = dis.T
        uniBest = ChromBest[np.where(dis != 0)[0]] # 相邻距离不为0的即为互异不重复的最优解
    else:
        uniBest = ChromBest
    
    return uniBest

def mintemp1(AIM_M, AIM_F, PUN_M, PUN_F, FieldDR, problem, maxormin, MAXGEN, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, distribute, drawing = 1):
    
    """
语法：
    该函数除了参数drawing外，不设置可缺省参数。当某个参数需要缺省时，在调用函数时传入None即可。
    比如当没有罚函数时，则在调用编程模板时将第3、4个参数设置为None即可，如：
    mintemp1(AIM_M, 'aimfuc', None, None, ..., maxormin)

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
    
    NIND : int - 种群规模，即种群中包含多少个个体
    
    SUBPOP : int - 子种群数量，即对一个种群划分多少个子种群
    
    GGAP : float - 代沟，表示子代与父代染色体及性状不相同的概率
    
    selectStyle : str - 指代所采用的低级选择算子的名称，如'rws'(轮盘赌选择算子)
    
    recombinStyle: str - 指代所采用的低级重组算子的名称，如'xovsp'(单点交叉)
    
    recopt : float - 交叉概率
    
    pm : float - 重组概率
    
    distribute : bool - 是否增强种群的分布性（可能会造成收敛慢）
    
    drawing : int - (可选参数)，0表示不绘图，1表示绘制最终结果图。默认drawing为1

输出参数：
    pop_trace : array - 种群进化记录器(进化追踪器),
                        第0列记录着各代种群最优个体的目标函数值
                        第1列记录着各代种群的适应度均值
                        第2列记录着各代种群最优个体的适应度值
    
    var_trace : array - 变量记录器，记录着各代种群最优个体的变量值，每一列对应一个控制变量
    
    times     : float - 进化所用时间

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
    """==========================初始化配置==========================="""
    # 获取目标函数和罚函数
    aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
    if PUN_F is not None:
        punishing = getattr(PUN_M, PUN_F) # 获得罚函数
    NVAR = FieldDR.shape[1] # 得到控制变量的个数
    # 定义进化记录器，初始值为nan
    pop_trace = (np.zeros((MAXGEN ,2)) * np.nan)
    # 定义变量记录器，记录方程解，初始长度为0
    sol_trace = np.zeros((0 ,NVAR))
    repnum = 0 # 初始化重复个体数为0
    """=========================开始遗传算法进化======================="""
    if problem == 'R':
        Chrom = ga.crtrp(NIND, FieldDR) # 生成初始种群
    elif problem == 'I':
        Chrom = ga.crtip(NIND, FieldDR)
    LegV = np.ones((NIND, 1)) # 初始化种群的可行性列向量
    [ObjV, LegV] = aimfuc(Chrom, LegV) # 求种群的目标函数值
    gen = 0
    badCounter = 0 # 用于记录在“遗忘策略下”被忽略的代数
    # 开始进化！！
    start_time = time.time() # 开始计时
    while gen < MAXGEN:
        if badCounter >= 10 * MAXGEN: # 若多花了10倍的迭代次数仍没有可行解出现，则跳出
            break
        FitnV = ga.ranking(maxormin * ObjV, LegV, None, SUBPOP)
        if PUN_F is not None:
            FitnV = punishing(LegV, FitnV) # 调用罚函数
        # 记录进化过程
        bestIdx = np.argmax(FitnV) # 获取最优个体的下标
        if LegV[bestIdx] != 0:
            feasible = np.where(LegV != 0)[0] # 排除非可行解
            pop_trace[gen,0] = np.sum(ObjV[feasible]) / ObjV[feasible].shape[0] # 记录种群个体平均目标函数值
            pop_trace[gen,1] = ObjV[bestIdx] # 记录当代目标函数的最优值
            allBestIdx = np.where(ObjV == ObjV[bestIdx])[0] # 找到相同的最优解
            # 获得互异不重复的最优个体
            ChromBest = Chrom[allBestIdx, :] # 取所有最优解，得到当代互异不重复的最优解
            uniBest = getUniBest(ChromBest) # 去除所有最优解中重复的解
            repnum = ChromBest.shape[0] - uniBest.shape[0] # 计算最优个体重复数
            sol_trace = np.vstack([sol_trace, uniBest]) # 把当代互异不重复最优解加入到sol_trace中
            badCounter = 0 # badCounter计数器清零
        else:
            gen -= 1 # 忽略这一代（遗忘策略）
            badCounter += 1
        if distribute == True: # 若要增强种群的分布性（可能会造成收敛慢）
            idx = np.argsort(ObjV[:, 0], 0)
            dis = np.diff(ObjV[idx,0]) / (np.max(ObjV[idx,0]) - np.min(ObjV[idx,0]) + 1)# 差分计算距离的修正偏移量
            dis = np.hstack([dis, dis[-1]])
            dis = dis + np.min(dis) # 修正偏移量+最小量=修正绝对量
            FitnV[idx, 0] *= np.exp(dis) # 根据相邻距离修改适应度，突出相邻距离大的个体，以增加种群的多样性
        # 进行遗传算子
        SelCh=ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP) # 选择
        SelCh=ga.recombin(recombinStyle, SelCh, recopt, SUBPOP) # 对所选个体进行重组
        if problem == 'R':
            SelCh=ga.mutbga(SelCh,FieldDR, pm) # 变异
            if repnum > Chrom.shape[0] * 0.01: # 当最优个体重复率高达1%时，进行一次高斯变异
                SelCh=ga.mutgau(SelCh, FieldDR, pm) # 高斯变异
        elif problem == 'I':
            SelCh=ga.mutint(SelCh, FieldDR, pm)
        LegVSel = np.ones((SelCh.shape[0], 1)) # 初始化育种种群的可行性列向量
        [ObjVSel, LegVSel] = aimfuc(SelCh, LegVSel) # 求育种种群的目标函数值
        FitnVSel = ga.ranking(maxormin * ObjVSel, LegVSel, None, SUBPOP) # 计算育种种群的适应度
        if PUN_F is not None:
            FitnVSel = punishing(LegVSel, FitnVSel) # 调用罚函数
        # 重插入
        [Chrom, ObjV, LegV] = ga.reins(Chrom, SelCh, SUBPOP, 1, 1, FitnV, FitnVSel, ObjV, ObjVSel, LegV, LegVSel)
        gen += 1
    end_time = time.time() # 结束计时
    times = end_time - start_time
    # 后处理进化记录器
    delIdx = np.where(np.isnan(pop_trace))[0]
    pop_trace = np.delete(pop_trace, delIdx, 0)
    if pop_trace.shape[0] == 0:
        raise RuntimeError('error: no feasible solution. (有效进化代数为0，没找到可行解。)')
    # 获得互异不重复的解
    sol_trace = getUniBest(sol_trace) # 去除所有最优解中重复的解
    # 此时sol_trace实际上会存在许多“假的”解，即很接近解，但是不是真正的解，因此需要后续处理
    tmp = np.ones((sol_trace.shape[0], 1)) # 无用参数，仅为凑够目标函数定义的输入参数
    sol_ObjVs = aimfuc(sol_trace, tmp)[0] # 计算解集对应的目标函数值
    best_solutions = sol_trace[[np.argmin(sol_ObjVs)], :] # 获得所有解中最优的即最接近真实解的一个解
    best_solutions_ObjV = aimfuc(best_solutions, np.ones((1, 1)))[0] # 计算最优的一个解对应的目标函数值
    sameBestIdx = np.where(sol_ObjVs == best_solutions_ObjV)[0]
    solutions = sol_trace[sameBestIdx, :] # 获得所有与最优解一样的解
    solutions_ObjVs = sol_ObjVs[sameBestIdx, :] # 获得所有与最优解一样的解对应的目标函数值
    # 绘图
    if drawing != 0:
        ga.trcplot(pop_trace, [['种群个体平均目标函数值', '种群最优个体目标函数值']])
    # 输出结果
    if solutions.shape[0] > 0:
        if best_solutions_ObjV == 0:
            print('方程的解为：')
        else: # 若best_solutions_ObjV不为0，则为近似解
            print('方程的近似解为：')
        for irun in range(solutions.shape[0]):
            print('解%s: '%(irun + 1))
            for i in range(NVAR):
                print('变量%d: %s'%((i + 1), solutions[irun][i]))
            print('对应的目标函数值为: %s'%(solutions_ObjVs[irun][0]))
            print('--------')
    print('有效进化代数：%s'%(pop_trace.shape[0]))
    print('时间已过 %s 秒'%(times))
    # 返回进化记录器、变量记录器以及执行时间
    return [pop_trace, solutions, times]

# -*- coding: utf-8 -*-

import numpy as np
import geatpy as ga
import time

def sga_mpc_real_templet(AIM_M, AIM_F, PUN_M, PUN_F, FieldDRs, problem, maxormin, MAXGEN, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, distribute, drawing = 1):
    
    """
sga_mpc_real_templet.py - 基于多种群竞争进化单目标编程模板(实值编码)

基于多种群竞争进化单目标编程模板(实值编码)，多种群之间加入竞争，同时各种群独立进化

语法：
    该函数除了drawing外，不设置可缺省参数。当某个参数需要缺省时，在调用函数时传入None即可。
    比如当没有罚函数时，则在调用编程模板时将第3、4个参数设置为None即可，如：
    sga_mpc_real_templet(AIM_M, 'aimfuc', None, None, ..., maxormin)

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
    
    GGAP : float - 代沟，本模板中该参数为无用参数，仅为了兼容同类的其他模板而设
    
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
    GGAP = 0.5 # 因为父子合并后选择，因此要将代沟设为0.5以维持种群规模
    # 获取目标函数和罚函数
    aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
    if PUN_F is not None:
        punishing = getattr(PUN_M, PUN_F) # 获得罚函数
    NVAR = FieldDRs[0].shape[1] # 得到控制变量的个数
    # 定义全局进化记录器，初始值为nan
    pop_trace = (np.zeros((MAXGEN ,2)) * np.nan)
    pop_trace[:, 0] = 0
    # 定义变量记录器，记录控制变量值，初始值为nan
    var_trace = (np.zeros((MAXGEN ,NVAR)) * np.nan)
    ax = None # 存储上一帧图形
    """=========================开始遗传算法进化======================="""
    Chroms = [] # 存储所有子种群的染色体
    LegVs = [] # 存储所有子种群的可行性列向量
    ObjVs = [] # 存储所有子种群的目标函数值
    repnums = [0] * len(FieldDRs) # 初始化重复个体数为0
    for index in range(len(FieldDRs)):
        if problem == 'R':
            Chrom = ga.crtrp(NIND, FieldDRs[index]) # 生成初始种群
        elif problem == 'I':
            Chrom = ga.crtip(NIND, FieldDRs[index])
        # 初始化相关变量
        Chroms.append(Chrom)
        LegV = np.ones((NIND, 1)) # 初始化初代种群的可行性列向量
        [ObjV, LegV] = aimfuc(Chrom, LegV) # 求初始种群的目标函数值
        ObjVs.append(ObjV)
        LegVs.append(LegV)
    gen = 0
    badCounter = 0 # 用于记录在“遗忘策略下”被忽略的代数
    # 开始进化！！
    start_time = time.time() # 开始计时
    while gen < MAXGEN:
        if badCounter >= 10 * MAXGEN: # 若多花了10倍的迭代次数仍没有可行解出现，则跳出
            break
        allChroms = np.zeros((NVAR, 0)).T # 存储所有子种群父子合并的染色体
        allObjVs = np.zeros((1, 0)).T # 存储所有子种群父子合并的目标函数值
        allLegVs = np.zeros((1, 0)).T # 存储所有子种群父子合并的可行性列向量
        if problem == 'R':
            for index in range(len(FieldDRs)): # 遍历各个子种群
                Chrom = Chroms[index] # 取某个子种群的Chrom
                ObjV = ObjVs[index] # 取某个子种群的ObjV
                LegV = LegVs[index] # 取某个子种群的LegV
                # 进行遗传算子，生成子代
                SelCh = ga.recombin(recombinStyle, Chrom, recopt, SUBPOP) # 重组
                SelCh = ga.mutbga(SelCh,FieldDRs[index], pm) # 变异
                if repnums[index] > Chrom.shape[0] * 0.01: # 当最优个体重复率高达1%时，进行一次高斯变异
                    SelCh = ga.mutgau(SelCh, FieldDRs[index], pm) # 高斯变异
                LegVSel = np.ones((SelCh.shape[0], 1)) # 初始化育种种群的可行性列向量
                [ObjVSel, LegVSel] = aimfuc(SelCh, LegVSel)
                # 父子合并
                Chrom = np.vstack([Chrom, SelCh]) 
                ObjV = np.vstack([ObjV, ObjVSel])
                LegV = np.vstack([LegV, LegVSel])
                allChroms = np.vstack([allChroms, Chrom])
                allObjVs = np.vstack([allObjVs, ObjV])
                allLegVs = np.vstack([allLegVs, LegV])
        elif problem == 'I':
            for index in range(len(FieldDRs)):
                Chrom = Chroms[index]
                # 进行遗传算子，生成子代
                SelCh = ga.recombin(recombinStyle, Chrom, recopt, SUBPOP) # 重组
                SelCh = ga.mutint(SelCh,FieldDRs[index], pm) # 变异
                LegVSel = np.ones((SelCh.shape[0], 1)) # 初始化育种种群的可行性列向量
                [ObjVSel, LegVSel] = aimfuc(SelCh, LegVSel) # 计算育种种群的目标函数值
                # 父子合并
                Chrom = np.vstack([Chrom, SelCh]) 
                ObjV = np.vstack([ObjV, ObjVSel])
                LegV = np.vstack([LegV, LegVSel])
                allChroms = np.vstack([allChroms, Chrom])
                allObjVs = np.vstack([allObjVs, ObjV])
                allLegVs = np.vstack([allLegVs, LegV])
        allFitnVs = ga.ranking(maxormin * allObjVs, allLegVs, None, SUBPOP) # 计算所有子种群的适应度
        if PUN_F is not None:
            allFitnVs = punishing(allLegVs, allFitnVs) # 调用罚函数
        # 记录进化过程
        bestIdx = np.argmax(allFitnVs)
        if allLegVs[bestIdx] != 0:
            feasible = np.where(allLegVs != 0)[0] # 排除非可行解
            pop_trace[gen,0] = np.sum(allObjVs[feasible]) / allObjVs[feasible].shape[0]
            pop_trace[gen,1] = allObjVs[bestIdx]
            var_trace[gen,:] = allChroms[bestIdx]
            # 绘制动态图
            if drawing == 2:
                ax = ga.sgaplot(pop_trace[:,[1]],'种群最优个体目标函数值', False, ax, gen)
            badCounter = 0 # badCounter计数器清零
        else:
            gen -= 1 # 忽略这一代（遗忘策略）
            badCounter += 1
        # 最后对合并的种群进行适应度评价并选出一半个体留到下一代（加入竞争的适应度评价）
        for index in range(len(FieldDRs)):
            Chrom = allChroms[index * NIND : (index + 2) * NIND]
            FitnV = allFitnVs[index * NIND : (index + 2) * NIND]
            ObjV = allObjVs[index * NIND : (index + 2) * NIND]
            LegV = allLegVs[index * NIND : (index + 2) * NIND]
            repnums[index] = len(np.where(FitnV[np.argmax(FitnV)] == FitnV)[0]) # 计算最优个体重复数
            if distribute == True: # 若要增强种群的分布性（可能会造成收敛慢）
                idx = np.argsort(ObjV[:, 0], 0)
                dis = np.diff(ObjV[idx,0]) / (np.max(ObjV[idx,0]) - np.min(ObjV[idx,0]) + 1)# 差分计算距离的修正偏移量
                dis = np.hstack([dis, dis[-1]])
                dis = dis + np.min(dis) # 修正偏移量+最小量=修正绝对量
                FitnV[idx, 0] *= np.exp(dis) # 根据相邻距离修改适应度，突出相邻距离大的个体，以增加种群的多样性
            [Chroms[index], ObjVs[index], LegVs[index]] = ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP, ObjV, LegV) # 选择
        gen += 1
    end_time = time.time() # 结束计时
    times = end_time - start_time
    # 后处理进化记录器
    delIdx = np.where(np.isnan(pop_trace))[0]
    pop_trace = np.delete(pop_trace, delIdx, 0)
    var_trace = np.delete(var_trace, delIdx, 0)
    if pop_trace.shape[0] == 0:
        raise RuntimeError('error: no feasible solution. (有效进化代数为0，没找到可行解。)')
    # 绘图
    if drawing != 0:
        ga.trcplot(pop_trace, [['种群个体平均目标函数值', '种群最优个体目标函数值']])
    # 输出结果
    if maxormin == 1:
        best_gen = np.argmin(pop_trace[:, 1]) # 记录最优种群是在哪一代
        best_ObjV = np.min(pop_trace[:, 1])
    elif maxormin == -1:
        best_gen = np.argmax(pop_trace[:, 1]) # 记录最优种群是在哪一代
        best_ObjV = np.max(pop_trace[:, 1])
    print('最优的目标函数值为：%s'%(best_ObjV))
    print('最优的控制变量值为：')
    for i in range(NVAR):
        print(var_trace[best_gen, i])
    print('有效进化代数：%s'%(pop_trace.shape[0]))
    print('最优的一代是第 %s 代'%(best_gen + 1))
    print('时间已过 %s 秒'%(times))
    # 返回进化记录器、变量记录器以及执行时间
    return [pop_trace, var_trace, times]

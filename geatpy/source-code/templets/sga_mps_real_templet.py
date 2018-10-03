# -*- coding: utf-8 -*-

import numpy as np
import geatpy as ga
import time

def sga_mps_real_templet(AIM_M, AIM_F, PUN_M, PUN_F, FieldDRs, problem, maxormin, MAXGEN, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, drawing = 1):
    
    """
sga_mps_real_templet.py - 基于多种群独立进化单目标编程模板(实值编码)

基于多种群独立进化单目标编程模板(实值编码)，各种群独立将父子两代合并进行选择，采取精英保留机制

语法：
    该函数除了drawing外，不设置可缺省参数。当某个参数需要缺省时，在调用函数时传入None即可。
    比如当没有罚函数时，则在调用编程模板时将第3、4个参数设置为None即可，如：
    sga_mps_real_templet(AIM_M, 'aimfuc', None, None, ..., maxormin)

输入参数：
    AIM_M - 目标函数的地址，由AIM_M = __import__('目标函数所在文件名')语句得到
            目标函数规范定义：f = aimfuc(Phen)
            其中Phen是种群的表现型矩阵
    
    AIM_F : str - 目标函数名
    
    PUN_M - 罚函数的地址，由PUN_M = __import__('罚函数所在文件名')语句得到
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
    
    NIND : int - 种群规模，即种群中包含多少个个体
    
    SUBPOP : int - 子种群数量，即对一个种群划分多少个子种群
    
    GGAP : float - 代沟，本模板中该参数为无用参数，仅为了兼容同类的其他模板而设
    
    selectStyle : str - 指代所采用的低级选择算子的名称，如'rws'(轮盘赌选择算子)
    
    recombinStyle: str - 指代所采用的低级重组算子的名称，如'xovsp'(单点交叉)
    
    recopt : float - 交叉概率
    
    pm : float - 重组概率
    
    drawing : int - (可选参数)，0表示不绘图，1表示绘制最终结果图。默认drawing为1

输出参数：
    pop_trace : array - 种群进化记录器(进化追踪器),
                        第0列记录着各代种群最优个体的目标函数值
                        第1列记录着各代种群的适应度均值
                        第2列记录着各代种群最优个体的适应度值
    
    var_trace : array - 变量记录器，记录着各代种群最优个体的变量值，每一列对应一个控制变量
    
    times     : float - 进化所用时间

模板使用注意：
    1.本模板调用的目标函数形如：aimfuc(Phen), 其中Phen表示种群的表现型矩阵
    2.本模板调用的罚函数形如: punishing(Phen, FitnV), 其中FitnV为用其他算法求得的适应度
      在罚函数定义中，必须将不满足约束条件的个体对应的适应度设为0，否则请修改模板使用

"""
 
    #==========================初始化配置==========================="""
    GGAP = 0.5 # 因为父子合并后选择，因此要将代沟设为0.5以维持种群规模
    # 获取目标函数和罚函数
    aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
    if PUN_F is not None:
        punishing = getattr(PUN_M, PUN_F) # 获得罚函数
    exIdx = np.array([]) # 存储非可行解的下标
    NVAR = FieldDRs[0].shape[1] # 得到控制变量的个数
    # 定义全局进化记录器，初始值为nan
    pop_trace = (np.zeros((MAXGEN ,2)) * np.nan)
    pop_trace[:, 0] = 0
    # 定义变量记录器，记录控制变量值，初始值为nan
    var_trace = (np.zeros((MAXGEN ,NVAR)) * np.nan) 
    """=========================开始遗传算法进化======================="""
    start_time = time.time() # 开始计时
    # 对于各个网格分别进行进化，采用全局进化记录器记录最优值
    for index in range(len(FieldDRs)):
        FieldDR = FieldDRs[index]
        if problem == 'R':
            Chrom = ga.crtrp(NIND, FieldDR) # 生成初始种群
        elif problem == 'I':
            Chrom = ga.crtip(NIND, FieldDR)
        repnum = 0 # 初始化重复个体数为0
        # 开始进化！！
        for gen in range(MAXGEN):
            # 进行遗传算子，生成子代
            SelCh=ga.recombin(recombinStyle, Chrom, recopt, SUBPOP) # 重组
            if problem == 'R':
                SelCh=ga.mutbga(SelCh,FieldDR, pm) # 变异
                if repnum > Chrom.shape[0] * 0.01: # 当最优个体重复率高达1%时，进行一次高斯变异
                    SelCh=ga.mutgau(SelCh, FieldDR, pm) # 高斯变异
            elif problem == 'I':
                SelCh=ga.mutint(SelCh, FieldDR, pm)
            Chrom = np.vstack([Chrom, SelCh]) # 父子合并
            ObjV = aimfuc(Chrom) # 求后代的目标函数值
            FitnV = ga.ranking(maxormin * ObjV, None, SUBPOP) # 适应度评价
            if PUN_F is not None:
                [FitnV, exIdx] = punishing(Chrom, FitnV) # 调用罚函数
            repnum = len(np.where(ObjV[np.argmax(FitnV)] == ObjV)[0]) # 计算最优个体重复数
            # 记录进化过程
            bestIdx = np.argmax(FitnV)
            wrongSign = np.ones((FitnV.shape[0], 1))
            wrongSign[list(exIdx)] = 0 # 对非可行解作标记
            if FitnV[bestIdx] != 0:
                if (np.isnan(pop_trace[gen,1])) or ((maxormin == 1) & (pop_trace[gen,1] >= ObjV[bestIdx])) or ((maxormin == -1) & (pop_trace[gen,1] <= ObjV[bestIdx])):
                    feasible = np.where(wrongSign != 0)[0] # 排除非可行解
                    pop_trace[gen,0] += np.sum(ObjV[feasible]) / ObjV[feasible].shape[0] / len(FieldDRs) # 记录种群个体平均目标函数值
                    pop_trace[gen,1] = ObjV[bestIdx] # 记录当代目标函数的最优值
                    var_trace[gen,:] = Chrom[bestIdx, :] # 记录当代最优的控制变量值
            Chrom=ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP) # 选择
    end_time = time.time() # 结束计时
    # 后处理进化记录器
    delIdx = np.where(np.isnan(pop_trace))[0]
    pop_trace = np.delete(pop_trace, delIdx, 0)
    var_trace = np.delete(var_trace, delIdx, 0)
    # 绘图
    if drawing == 1:
        ga.trcplot(pop_trace, [['种群个体平均目标函数值', '种群最优个体目标函数值']])
    # 输出结果
    if maxormin == 1:
        best_gen = np.argmin(pop_trace[:, 1]) # 记录最优种群是在哪一代
        print('最优的目标函数值为：', np.min(pop_trace[:, 1]))
    elif maxormin == -1:
        best_gen = np.argmax(pop_trace[:, 1]) # 记录最优种群是在哪一代
        print('最优的目标函数值为：', np.max(pop_trace[:, 1]))
    print('最优的控制变量值为：')
    for i in range(NVAR):
        print(var_trace[best_gen, i])
    print('最优的一代是第', best_gen + 1, '代')
    times = end_time - start_time
    print('时间已过', times, '秒')
    # 返回进化记录器、变量记录器以及执行时间
    return [pop_trace, var_trace, times]
    
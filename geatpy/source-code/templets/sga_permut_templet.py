# -*- coding: utf-8 -*-

import numpy as np
import geatpy as ga
import time

def sga_permut_templet(AIM_M, AIM_F, PUN_M, PUN_F, NVAR, VarLen, maxormin, MAXGEN, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, distribute, drawing = 1):
    
    """
sga_permut_templet.py - 单目标编程模板(排列编码)

排列编码即每条染色体的基因都是无重复正整数的编码方式。

语法：
    该函数除了参数drawing外，不设置可缺省参数。当某个参数需要缺省时，在调用函数时传入None即可。
    比如当没有罚函数时，则在调用编程模板时将第3、4个参数设置为None即可，如：
    sga_permut_templet(AIM_M, 'aimfuc', None, None, ..., maxormin)

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
    
    NVAR : int - 变量个数，排列编码的染色体长度等于变量个数
    
    VarLen : int - 排列集合的大小
                   例如VarLen = 5表示是从1,2,3,4,5中抽取若干个数排列组成染色体
    
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
    # 定义进化记录器，初始值为nan
    pop_trace = (np.zeros((MAXGEN ,2)) * np.nan)
    # 定义变量记录器，记录控制变量值，初始值为nan
    var_trace = (np.zeros((MAXGEN ,NVAR)) * np.nan) 
    ax = None # 存储上一帧图形
    """=========================开始遗传算法进化======================="""
    #生成初始种群
    Chrom = ga.crtpp(NIND, NVAR, VarLen)
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
            var_trace[gen,:] = Chrom[bestIdx, :] # 记录当代最优的控制变量值
            # 绘制动态图
            if drawing == 2:
                ax = ga.sgaplot(pop_trace[:,[1]],'种群最优个体目标函数值', False, ax, gen)
            badCounter = 0 # badCounter计数器清零
        else:
            gen -= 1 # 忽略这一代（遗忘策略）
            badCounter += 1
        # 进行遗传算子
        SelCh=ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP) # 选择
        SelCh=ga.recombin(recombinStyle, Selch, recopt, SUBPOP) # 对选择的个体进行重组
        SelCh=ga.mutpp(SelCh, VarLen, pm) # 排列编码种群变异
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
    print('最优的目标函数值为：%f'%(best_ObjV))
    print('最优的控制变量值为：')
    for i in range(NVAR):
        print(var_trace[best_gen, i])
    print('有效进化代数：%d'%(pop_trace.shape[0]))
    print('最优的一代是第 %d 代'%(best_gen + 1))
    print('时间已过 %f 秒'%(times))
    # 返回进化记录器、变量记录器以及执行时间
    return [pop_trace, var_trace, times]

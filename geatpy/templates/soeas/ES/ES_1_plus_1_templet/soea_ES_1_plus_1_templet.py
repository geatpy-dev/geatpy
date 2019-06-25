# -*- coding: utf-8 -*-
import time
import numpy as np
import geatpy as ea # 导入geatpy库
from sys import path as paths
from os import path
paths.append(path.split(path.split(path.realpath(__file__))[0])[0])

class soea_ES_1_plus_1_templet(ea.Algorithm):
    
    """
soea_ES_1_plus_1_templet : class - (1+1)进化策略模板

算法描述:
    本模板实现的是(1+1)进化策略。算法流程如下：
    1) 根据编码规则初始化N个个体的种群。
    2) 若满足停止条件则停止，否则继续执行。
    3) 对当前种群进行统计分析，比如记录其最优个体、平均适应度等等。
    4) 初始化变异压缩概率Shrink，该压缩概率用于控制高斯变异中的标准差。
    5) 独立地对这种群个体进行高斯变异，得到试验种群。
    6) 在当前种群和实验种群之间采用一对一生存者选择方法得到新一代种群，
    同时统计新一代种群中有多少个个体继承自实验种群（即变异成功率）。
    7) 根据变异成功率修改压缩概率Shrink。
    8) 回到第2步。
    为了使算法也能很好地求解约束优化问题，本算法模板稍作修改，增添“遗忘策略”，
    当某一代没有可行个体时，让进化记录器忽略这一代，不对这一代的个体进行记录。
    
模板使用注意：
    1.本模板调用的目标函数形如：[ObjV,CV] = aimFuc(Vars, CV), 
      其中Vars表示决策变量矩阵, CV为种群的违反约束程度矩阵(详见Geatpy数据结构),
    若不符合上述规范，则请修改算法模板或自定义新算法模板。

参考文献:
    [1] Beyer H G , Schwefel H P . Evolution strategies – A comprehensive 
    introduction[J]. Natural Computing, 2002, 1(1):3-52.

"""
    
    def __init__(self, problem, population):
        self.name = '(1+1)ES'
        self.problem = problem
        self.population = population
        if population.Encoding == 'R' or population.Encoding == 'I':
            self.mutFunc = 'mutgau' # 高斯变异
        self.drawing = 1 # 绘图
        self.ax = None # 存储上一桢动画
        self.passTime = 0 # 记录用时
    
    def stat(self, population): # 分析记录
        # 进行进化记录
        feasible = np.where(np.all(population.CV <= 0, 1))[0] # 找到可行解的下标
        if len(feasible) > 0:
            tempPop = population[feasible]
            bestIdx = np.argmax(tempPop.FitnV) # 获取最优个体的下标 
            self.obj_trace[self.currentGen,0] = np.sum(tempPop.ObjV) / tempPop.sizes # 记录种群个体平均目标函数值
            self.obj_trace[self.currentGen,1] = tempPop.ObjV[bestIdx] # 记录当代目标函数的最优值
            self.var_trace[self.currentGen,:] = tempPop.Phen[bestIdx, :] # 记录当代最优的决策变量值
            self.forgetCount = 0 # “遗忘策略”计数器清零
            self.passTime += time.time() - self.timeSlot # 更新用时记录
            if self.drawing == 2:
                self.ax = ea.soeaplot(self.obj_trace[:,[1]],'种群最优个体目标函数值', False, self.ax, self.currentGen) # 绘制动态图
            self.timeSlot = time.time() # 更新时间戳
        else:
            self.currentGen -= 1 # 忽略这一代
            self.forgetCount += 1 # “遗忘策略”计数器加1
    
    def terminated(self, population): # 判断是否终止进化
        if self.currentGen < self.MAXGEN or self.forgetCount >= 10 * self.MAXGEN:
            self.stat(population) # 分析记录当代种群的数据
            self.currentGen += 1 # 进化代数+1
            return False
        else:
            return True
    
    def run(self):
        #==========================初始化配置===========================
        population = self.population
        NIND = population.sizes
        NVAR = self.problem.Dim # 得到决策变量的个数
        self.obj_trace = (np.zeros((self.MAXGEN, 2)) * np.nan) # 定义目标函数值记录器，初始值为nan
        self.var_trace = (np.zeros((self.MAXGEN, NVAR)) * np.nan) # 定义变量记录器，记录决策变量值，初始值为nan
        self.forgetCount = 0 # “遗忘策略”计数器，用于记录连续出现最优个体不是可行个体的代数
        #===========================准备进化============================
        self.timeSlot = time.time() # 开始计时
        if population.Chrom is None:
            population.initChrom(NIND) # 初始化种群染色体矩阵（内含染色体解码，详见Population类的源码）
        population.ObjV, population.CV = self.problem.aimFuc(population.Phen, population.CV) # 计算种群的目标函数值
        population.FitnV = ea.scaling(self.problem.maxormins * population.ObjV, population.CV) # 计算适应度
        self.evalsNum = population.sizes # 记录评价次数
        Sigma = 0.5 * (population.Field[1,:] - population.Field[0,:]) / 3 # 初始化高斯变异的Sigma
        #===========================开始进化============================
        self.currentGen = 0
        while self.terminated(population) == False:
            # 进行进化操作
            experimentPop = population.copy() # 存储试验种群
            experimentPop.Chrom = ea.mutate('mutgau', experimentPop.Encoding, experimentPop.Chrom, experimentPop.Field, experimentPop.Lind, Sigma) # 变异（这里变异概率设为染色体长度）
            # 求进化后个体的目标函数值
            experimentPop.Phen = experimentPop.decoding() # 染色体解码
            experimentPop.ObjV, experimentPop.CV = self.problem.aimFuc(experimentPop.Phen, experimentPop.CV)
            self.evalsNum += population.sizes # 更新评价次数
            tempPop = population + experimentPop # 临时合并，以调用otos进行一对一生存者选择
            tempPop.FitnV = ea.scaling(self.problem.maxormins * tempPop.ObjV, tempPop.CV) # 计算适应度
            chooseIdx = ea.selecting('otos', tempPop.FitnV, NIND) # 采用One-to-One Survivor选择
            population = tempPop[chooseIdx] # 产生新一代种群
            # 利用1/5规则调整变异压缩概率（实质上是通过变异压缩概率来调整高斯变异的标准差，详见mutgau帮助文档）
            successfulRate = len(np.where(chooseIdx >= NIND)[0]) / (2 * NIND)
            if successfulRate < 1/5:
                Sigma *= 0.817
            elif successfulRate > 1/5:
                Sigma /= 0.817
        # 处理进化记录器
        delIdx = np.where(np.isnan(self.obj_trace))[0]
        self.obj_trace = np.delete(self.obj_trace, delIdx, 0)
        self.var_trace = np.delete(self.var_trace, delIdx, 0)
        if self.obj_trace.shape[0] == 0:
            raise RuntimeError('error: No feasible solution. (有效进化代数为0，没找到可行解。)')
        self.passTime += time.time() - self.timeSlot # 更新用时记录
        # 绘图
        if self.drawing != 0:
            ea.trcplot(self.obj_trace, [['种群个体平均目标函数值', '种群最优个体目标函数值']])
        # 返回最后一代种群、进化记录器、变量记录器以及执行时间
        return [population, self.obj_trace, self.var_trace]

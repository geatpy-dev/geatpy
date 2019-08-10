# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea # 导入geatpy库
from sys import path as paths
from os import path
paths.append(path.split(path.split(path.realpath(__file__))[0])[0])

class soea_ES_1_plus_1_templet(ea.SoeaAlgorithm):
    
    """
soea_ES_1_plus_1_templet : class - (1+1)进化策略模板

算法描述:
    本模板实现的是(1+1)进化策略。算法流程如下：
    1) 根据编码规则初始化N个个体的种群。
    2) 若满足停止条件则停止，否则继续执行。
    3) 对当前种群进行统计分析，比如记录其最优个体、平均适应度等等。
    4) 初始化控制高斯变异中的标准差Sigma。
    5) 独立地对这种群个体进行高斯变异，得到试验种群。
    6) 在当前种群和实验种群之间采用一对一生存者选择方法得到新一代种群，
    同时统计新一代种群中有多少个个体继承自实验种群（即变异成功率）。
    7) 根据变异成功率修改Sigma。
    8) 回到第2步。

模板使用注意:
    本模板调用的目标函数形如：aimFunc(pop), 
    其中pop为Population类的对象，代表一个种群，
    pop对象的Phen属性（即种群染色体的表现型）等价于种群所有个体的决策变量组成的矩阵，
    该函数根据该Phen计算得到种群所有个体的目标函数值组成的矩阵，并将其赋值给pop对象的ObjV属性。
    若有约束条件，则在计算违反约束程度矩阵CV后赋值给pop对象的CV属性（详见Geatpy数据结构）。
    该函数不返回任何的返回值，求得的目标函数值保存在种群对象的ObjV属性中，
                          违反约束程度矩阵保存在种群对象的CV属性中。
    例如：population为一个种群对象，则调用aimFunc(population)即可完成目标函数值的计算，
         此时可通过population.ObjV得到求得的目标函数值，population.CV得到违反约束程度矩阵。
    若不符合上述规范，则请修改算法模板或自定义新算法模板。

参考文献:
    [1] Beyer H G , Schwefel H P . Evolution strategies – A comprehensive 
    introduction[J]. Natural Computing, 2002, 1(1):3-52.

"""
    
    def __init__(self, problem, population):
        ea.SoeaAlgorithm.__init__(self, problem, population) # 先调用父类构造方法
        if str(type(population)) != "<class 'Population.Population'>":
            raise RuntimeError('传入的种群对象必须为Population类型')
        self.name = '(1+1)ES'
        if population.Encoding == 'RI':
            self.mutFunc = 'mutgau' # 高斯变异
        else:
            raise RuntimeError('编码方式必须为''RI''.')
            
    def run(self):
        #==========================初始化配置===========================
        population = self.population
        NIND = population.sizes
        self.initialization() # 初始化算法模板的一些动态参数
        #===========================准备进化============================
        if population.Chrom is None:
            population.initChrom(NIND) # 初始化种群染色体矩阵（内含染色体解码，详见Population类的源码）
        self.problem.aimFunc(population) # 计算种群的目标函数值
        population.FitnV = ea.scaling(self.problem.maxormins * population.ObjV, population.CV) # 计算适应度
        self.evalsNum = population.sizes # 记录评价次数
        Sigma = 0.5 * (population.Field[1,:] - population.Field[0,:]) / 3 # 初始化高斯变异的Sigma
        #===========================开始进化============================
        while self.terminated(population) == False:
            # 进行进化操作
            experimentPop = population.copy() # 存储试验种群
            experimentPop.Chrom = ea.mutate('mutgau', experimentPop.Encoding, experimentPop.Chrom, experimentPop.Field, experimentPop.Lind, Sigma) # 变异（这里变异概率设为染色体长度）
            # 求进化后个体的目标函数值
            experimentPop.Phen = experimentPop.decoding() # 染色体解码
            self.problem.aimFunc(experimentPop) # 计算目标函数值
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
        
        return self.finishing(population) # 调用finishing完成后续工作并返回结果

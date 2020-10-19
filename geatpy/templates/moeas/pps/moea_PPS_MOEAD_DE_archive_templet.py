# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea  # 导入geatpy库
from scipy.spatial.distance import cdist
from sys import path as paths
from os import path

paths.append(path.split(path.split(path.realpath(__file__))[0])[0])


class moea_PPS_MOEAD_DE_archive_templet(ea.MoeaAlgorithm):
    """
moea_PPS_MOEAD_DE_archive_templet : class - 基于pps策略的带全局存档的多目标进化MOEA/D-DE算法模板
    
算法描述:
    采用PPS-MOEA/D-DE进行多目标优化，PPS策略详见参考文献[1]，
    注：MOEA/D不适合在Python上实现，在Python上，MOEA/D的性能会大幅度降低。

参考文献:
    [1] Zhun Fan, Wenji Li, Xinye Cai*, Hui Li, Caimin Wei, Qingfu Zhang, 
    Kalyanmoy Deb, and Erik Goodman. Push and Pull Search for Solving 
    Constrained Multi-objective Optimization Problems, Swarm and Evolutionary 
    Computation, vol. 44, no. 2, pp. 665-679, 2019.
    
    """

    def __init__(self, problem, population):
        ea.MoeaAlgorithm.__init__(self, problem, population)  # 先调用父类构造方法
        if population.ChromNum != 1:
            raise RuntimeError('传入的种群对象必须是单染色体的种群类型。')
        self.name = 'PPS-MOEA/D-DE-archive'
        if population.Encoding == 'RI':
            self.F = 0.5  # DE的F
            self.Cr = 1.0  # DE的Cr
            self.mutOper = ea.Mutpolyn(Pm=1 / self.problem.Dim, DisI=20)  # 生成多项式变异算子对象
        else:
            raise RuntimeError('编码方式必须为''RI''.')
        self.neighborSize = None  # 邻域大小，当设置为None时，将会自动设置为等于种群规模的十分之一
        if self.problem.M <= 2:
            self.decomposition = ea.tcheby  # 采用切比雪夫权重聚合法
        else:
            self.decomposition = ea.pbi  # 采用pbi权重聚合法
        self.Ps = 0.9  # (Probability of Selection)表示进化时有多大的概率只从邻域中选择个体参与进化
        self.Nr = 2  # MOEAD-DE中的参数nr，默认为2
        self.MAXSIZE = population.sizes  # 全局非支配解存档的大小限制，这里设为等于初始设定的种群个体数
        # PPS策略的一些需要设置的参数
        self.Tc = 0.8  # 论文中的Tc，这里暂设为0.8，在run()函数中它将乘上MAXGEN
        self.LastLGen = 20  # 论文中的参数l
        self.varient_epsilon = 1e-3  # 论文中的参数varient_epsilon
        self.alpha = 0.95  # 论文中的α
        self.tao = 0.1  # 论文中的𝜏
        self.cp = 2  # 论文中的cp

    def create_offspring(self, population, Xr0, select_rand, Mask, neighbor_index, idealPoint):

        """
        描述:
            该函数用于产生子代个体以及更新理想点，它实际上是下面的主代码里抽取出来的，
            若有理解困难，可以把该函数的代码重新放入主代码中。
        """
        if select_rand < self.Ps:
            indices = neighbor_index
        else:
            indices = np.arange(population.sizes)
        offspring = ea.Population(population.Encoding, population.Field, 1)  # 实例化一个种群对象用于存储进化的后代（这里只进化生成一个后代）
        r = indices[ea.rps(len(indices), 2)]  # 随机选择两个索引作为差分向量的索引
        r1, r2 = r[0], r[1]  # 得到差分向量索引
        offspring.Chrom = Xr0
        offspring.Chrom[0][Mask] = offspring.Chrom[0][Mask] + self.F * (
                    population.Chrom[r1][Mask] - population.Chrom[r2][Mask])
        offspring.Chrom = self.mutOper.do(offspring.Encoding, offspring.Chrom, offspring.Field)  # 多项式变异
        self.call_aimFunc(offspring)  # 求进化后个体的目标函数值
        # 更新理想点
        idealPoint = ea.crtidp(offspring.ObjV, maxormins=self.problem.maxormins, old_idealPoint=idealPoint)
        return offspring, indices, idealPoint

    def push_stage_reinsertion(self, indices, population, offspring, idealPoint, referPoint):

        """
        描述:
            适用于push stage的重插入更新种群个体。
        """

        weights = referPoint[indices, :]
        pop_ObjV = population.ObjV[indices, :]  # 获取邻居个体的目标函数值
        CombinObjV = self.decomposition(pop_ObjV, weights, idealPoint, maxormins=self.problem.maxormins)
        off_CombinObjV = self.decomposition(offspring.ObjV, weights, idealPoint, maxormins=self.problem.maxormins)
        population[indices[np.where(off_CombinObjV <= CombinObjV)[0][:self.Nr]]] = offspring

    def pull_stage_reinsertion(self, indices, population, offspring, idealPoint, referPoint, epsilon_k):

        """
        描述:
            适用于pull stage的重插入更新种群个体。
        """

        weights = referPoint[indices, :]
        pop_ObjV = population.ObjV[indices, :]  # 获取邻居个体的目标函数值
        CombinObjV = self.decomposition(pop_ObjV, weights, idealPoint, maxormins=self.problem.maxormins)
        off_CombinObjV = self.decomposition(offspring.ObjV, weights, idealPoint, maxormins=self.problem.maxormins)
        Violation = ea.mergecv(population.CV[indices, :])
        off_Violation = ea.mergecv(offspring.CV)
        population[(indices[np.where((off_CombinObjV <= CombinObjV) &
                                     ((Violation <= epsilon_k) & (off_Violation <= epsilon_k) | (
                                                 Violation == off_Violation)) |
                                     (off_Violation < Violation))[0]])[:self.Nr]] = offspring

    def updateNDSet(self, population, globalNDSet=None):

        """
        描述:
            更新globalNDSet。
        """

        if globalNDSet is None:
            globalNDSet = population
        else:
            globalNDSet = population + globalNDSet  # 将population与全局归档集合并
        if globalNDSet.CV is not None:  # CV不为None说明有设置约束条件
            globalNDSet = globalNDSet[np.where(np.all(globalNDSet.CV <= 0, 1))[0]]  # 排除非可行解
        if globalNDSet.sizes != 0:
            [levels, criLevel] = ea.ndsortDED(globalNDSet.ObjV, None, None, globalNDSet.CV,
                                              self.problem.maxormins)  # 非支配排序
            globalNDSet = globalNDSet[np.where(levels == 1)[0]]
        if globalNDSet.sizes > self.MAXSIZE:
            dis = ea.crowdis(globalNDSet.ObjV, np.ones(globalNDSet.sizes))  # 计算拥挤距离
            globalNDSet = globalNDSet[np.argsort(-dis)[:self.MAXSIZE]]  # 根据拥挤距离选择符合个数限制的解保留在存档中
        return globalNDSet

    def run(self, prophetPop=None):  # prophetPop为先知种群（即包含先验知识的种群）
        # ==========================初始化配置===========================
        population = self.population
        self.initialization()  # 初始化算法模板的一些动态参数
        pushStage = True  # 一开始是push stage
        rk = 1.0  # 论文中的rk，k的含义在论文中是代数，这里保留名称不作变化，下同
        epsilon_k = 0  # 论文中的𝜀(k)
        epsilon_0 = 0  # 论文中的𝜀(0)
        idealPoints = []  # 存储历代的理想点的列表
        nadirPoints = []  # 存储历代的反理想点的列表
        delta = np.array([1e-6] * self.problem.M)  # 论文中为了避免分母为0而设的delta
        self.Tc *= self.MAXGEN
        self.LastLGen = min(self.LastLGen, self.MAXGEN)
        # ===========================准备进化============================
        uniformPoint, NIND = ea.crtup(self.problem.M, population.sizes)  # 生成在单位目标维度上均匀分布的参考点集
        population.initChrom(NIND)  # 初始化种群染色体矩阵，此时种群规模将调整为uniformPoint点集的大小，initChrom函数会把种群规模给重置
        self.call_aimFunc(population)  # 计算种群的目标函数值
        # 插入先验知识（注意：这里不会对先知种群prophetPop的合法性进行检查，故应确保prophetPop是一个种群类且拥有合法的Chrom、ObjV、Phen等属性）
        if prophetPop is not None:
            population = (prophetPop + population)[:NIND]  # 插入先知种群
        # 确定邻域大小
        if self.neighborSize is None:
            self.neighborSize = population.sizes // 10
        self.neighborSize = max(self.neighborSize, 2)  # 确保不小于2
        # 生成由所有邻居索引组成的矩阵
        neighborIdx = np.argsort(cdist(uniformPoint, uniformPoint), axis=1, kind='mergesort')[:, :self.neighborSize]
        # 计算理想点
        idealPoint = ea.crtidp(population.ObjV, maxormins=self.problem.maxormins)
        # 创建全局存档
        globalNDSet = self.updateNDSet(population)
        # ===========================开始进化============================
        while self.terminated(population) == False:
            idealPoints.append(idealPoint)
            nadirPoints.append(ea.crtidp(population.ObjV, maxormins=self.problem.maxormins, reverse=True))
            # 更新epsilon_k
            if self.currentGen < self.Tc:
                # 更新rk
                if self.currentGen >= self.LastLGen:
                    past_gen = self.currentGen - self.LastLGen
                    rk = np.max(
                        [np.abs((idealPoints[-1] - idealPoints[past_gen]) / np.max([idealPoints[past_gen], delta], 0)),
                         np.abs((nadirPoints[-1] - nadirPoints[past_gen]) / np.max([nadirPoints[past_gen], delta], 0))])
                violation, count = ea.mergecv(population.CV, return_count=True)
                if rk <= self.varient_epsilon and pushStage:
                    epsilon_0 = np.max(violation)
                    epsilon_k = epsilon_0
                    pushStage = False
                if pushStage == False:
                    rf = count / population.sizes
                    if rf < self.alpha:
                        epsilon_k *= (1 - self.tao)
                    else:
                        epsilon_k = (1 - self.currentGen / self.Tc) ** self.cp * epsilon_0
            else:
                epsilon_k = 0
            # 分开push stage和pull stage进行进化
            select_rands = np.random.rand(population.sizes)
            Masks = np.random.rand(population.sizes, population.Lind) < self.Cr
            if pushStage:
                for i in range(population.sizes):
                    # 产生后代
                    offspring, indices, idealPoint = self.create_offspring(population, population.Chrom[[i], :],
                                                                           select_rands[i], Masks[i], neighborIdx[i, :],
                                                                           idealPoint)
                    # 重插入
                    self.push_stage_reinsertion(indices, population, offspring, idealPoint, uniformPoint)  # 重插入更新种群个体
            else:
                for i in range(population.sizes):
                    # 产生后代
                    offspring, indices, idealPoint = self.create_offspring(population, population.Chrom[[i], :],
                                                                           select_rands[i], Masks[i], neighborIdx[i, :],
                                                                           idealPoint)
                    # 重插入
                    self.pull_stage_reinsertion(indices, population, offspring, idealPoint, uniformPoint, epsilon_k)
            # 完成当代的进化后，更新全局存档
            globalNDSet = self.updateNDSet(population, globalNDSet)
        return self.finishing(population, globalNDSet)  # 调用finishing完成后续工作并返回结果

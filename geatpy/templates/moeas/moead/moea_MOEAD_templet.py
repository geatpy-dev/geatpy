# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea  # 导入geatpy库
from scipy.spatial.distance import cdist
from sys import path as paths
from os import path

paths.append(path.split(path.split(path.realpath(__file__))[0])[0])


class moea_MOEAD_templet(ea.MoeaAlgorithm):
    """
moea_MOEAD_templet : class - 多目标进化MOEA/D算法模板（采用可行性法则处理约束）
    
算法描述:
    采用MOEA/D（不设全局最优存档）进行多目标优化，算法详见参考文献[1]。
    注：MOEA/D不适合在Python上实现，在Python上，MOEA/D的性能会大幅度降低。

参考文献:
    [1] Qingfu Zhang, Hui Li. MOEA/D: A Multiobjective Evolutionary Algorithm 
    Based on Decomposition[M]. IEEE Press, 2007.

    """

    def __init__(self, problem, population):
        ea.MoeaAlgorithm.__init__(self, problem, population)  # 先调用父类构造方法
        if population.ChromNum != 1:
            raise RuntimeError('传入的种群对象必须是单染色体的种群类型。')
        self.name = 'MOEA/D'
        if population.Encoding == 'P':
            self.recOper = ea.Xovpmx(XOVR=1, Half_N=True)  # 生成部分匹配交叉算子对象
            self.mutOper = ea.Mutinv(Pm=1)  # 生成逆转变异算子对象
        elif population.Encoding == 'BG':
            self.recOper = ea.Xovud(XOVR=1, Half_N=True)  # 生成均匀交叉算子对象
            self.mutOper = ea.Mutbin(Pm=None)  # 生成二进制变异算子对象，Pm设置为None时，具体数值取变异算子中Pm的默认值
        elif population.Encoding == 'RI':
            self.recOper = ea.Recsbx(XOVR=1, n=20, Half_N=True)  # 生成模拟二进制交叉算子对象
            self.mutOper = ea.Mutpolyn(Pm=1 / self.problem.Dim, DisI=20)  # 生成多项式变异算子对象
        else:
            raise RuntimeError('编码方式必须为''BG''、''RI''或''P''.')
        self.neighborSize = None  # 邻域大小，当设置为None时，将会自动设置为等于种群规模
        if self.problem.M <= 2:
            self.decomposition = ea.tcheby  # 采用切比雪夫权重聚合法
        else:
            self.decomposition = ea.pbi  # 采用pbi权重聚合法
        self.Ps = 0.9  # (Probability of Selection)表示进化时有多大的概率只从邻域中选择个体参与进化

    def reinsertion(self, indices, population, offspring, idealPoint, referPoint):

        """
        描述:
            重插入更新种群个体。
        """

        weights = referPoint[indices, :]
        pop_ObjV = population.ObjV[indices, :]  # 获取邻居个体的目标函数值
        pop_CV = population.CV[indices, :] if population.CV is not None else None  # 获取邻居个体的违反约束程度矩阵
        CombinObjV = self.decomposition(pop_ObjV, weights, idealPoint, pop_CV, self.problem.maxormins)
        off_CombinObjV = self.decomposition(offspring.ObjV, weights, idealPoint, offspring.CV, self.problem.maxormins)
        population[indices[np.where(off_CombinObjV <= CombinObjV)[0]]] = offspring

    def run(self, prophetPop=None):  # prophetPop为先知种群（即包含先验知识的种群）
        # ==========================初始化配置===========================
        population = self.population
        self.initialization()  # 初始化算法模板的一些动态参数
        # ===========================准备进化============================
        uniformPoint, NIND = ea.crtup(self.problem.M, population.sizes)  # 生成在单位目标维度上均匀分布的参考点集
        population.initChrom(NIND)  # 初始化种群染色体矩阵，此时种群规模将调整为uniformPoint点集的大小，initChrom函数会把种群规模给重置
        self.call_aimFunc(population)  # 计算种群的目标函数值
        # 插入先验知识（注意：这里不会对先知种群prophetPop的合法性进行检查，故应确保prophetPop是一个种群类且拥有合法的Chrom、ObjV、Phen等属性）
        if prophetPop is not None:
            population = (prophetPop + population)[:NIND]  # 插入先知种群
        # 确定邻域大小
        if self.neighborSize is None:
            self.neighborSize = population.sizes
        self.neighborSize = max(self.neighborSize, 2)  # 确保不小于2
        # 生成由所有邻居索引组成的矩阵
        neighborIdx = np.argsort(cdist(uniformPoint, uniformPoint), axis=1, kind='mergesort')[:, :self.neighborSize]
        # 计算理想点
        idealPoint = ea.crtidp(population.ObjV, population.CV, self.problem.maxormins)
        # ===========================开始进化============================
        while self.terminated(population) == False:
            select_rands = np.random.rand(population.sizes)  # 生成一组随机数
            for i in range(population.sizes):
                indices = neighborIdx[i, :]  # 得到邻居索引
                if select_rands[i] < self.Ps:
                    chooseIdx = indices[ea.rps(self.neighborSize, 2)]  # 只从邻域中选择
                else:
                    chooseIdx = ea.rps(population.sizes, 2)
                matting_Chrom = population.Chrom[chooseIdx, :]  # 选出2条来自被选个体的染色体
                offspring = ea.Population(population.Encoding, population.Field, 1)  # 实例化一个种群对象用于存储进化的后代（这里只进化生成一个后代）
                # 对选出的个体进行进化操作
                offspring.Chrom = self.recOper.do(matting_Chrom)  # 重组
                offspring.Chrom = self.mutOper.do(offspring.Encoding, offspring.Chrom, offspring.Field)  # 变异
                self.call_aimFunc(offspring)  # 求进化后个体的目标函数值
                # 更新理想点
                idealPoint = ea.crtidp(offspring.ObjV, offspring.CV, self.problem.maxormins, idealPoint)
                # 重插入更新种群个体
                self.reinsertion(indices, population, offspring, idealPoint, uniformPoint)
        return self.finishing(population)  # 调用finishing完成后续工作并返回结果

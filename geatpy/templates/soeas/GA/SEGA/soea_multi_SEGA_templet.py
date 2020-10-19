# -*- coding: utf-8 -*-
import geatpy as ea  # 导入geatpy库
import numpy as np
from sys import path as paths
from os import path

paths.append(path.split(path.split(path.realpath(__file__))[0])[0])


class soea_multi_SEGA_templet(ea.SoeaAlgorithm):
    """
soea_multi_SEGA_templet : class - Multi-population Strengthen Elitist GA templet(增强精英保留的多种群协同遗传算法模板)

模板说明:
    该模板是内置算法模板soea_SEGA_templet的多种群协同版本，为不同的种群设置不同的重组和变异概率。
    注意：本算法模板中的population为一个存储种群类对象的列表，而不是单个种群类对象。

算法描述:
    本模板实现的是增强精英保留的多种群协同遗传算法。算法流程如下：
    1) 循环population列表，初始化列表中的各个种群的染色体，并将所有种群所有个体的数目记录为NindAll。
    2) 若满足进化算法停止条件则停止，否则继续执行。
    3) 循环对各个种群独立进行选择、重组、变异，得到各个种群的子代，并将父代和子代个体合并。
    4) 对所有种群的所有个体进行统一的适应度评价。
    5) 根据适应度调用选择算子进行环境选择，选择出NindAll个个体形成新一代种群。
    6) 根据概率进行种群间个体迁移。
    7) 回到第2步。
    
"""

    def __init__(self, problem, population):
        ea.SoeaAlgorithm.__init__(self, problem, population)  # 先调用父类构造方法
        if type(population) != list:
            raise RuntimeError('传入的种群对象列表必须为list类型')
        self.name = 'multi-SEGA'
        self.PopNum = len(population)  # 种群数目
        self.selFunc = 'tour'  # 锦标赛选择算子
        self.migFr = 5  # 发生种群迁移的间隔代数
        self.migOpers = ea.Migrate(MIGR=0.2, Structure=2, Select=1, Replacement=2)  # 生成种群迁移算子对象
        # 为不同的种群设置不同的重组、变异算子
        self.recOpers = []
        self.mutOpers = []
        Pms = np.linspace(1 / self.problem.Dim, 1, self.PopNum)  # 生成变异概率列表，为不同的种群分配不同的变异概率
        Pcs = np.linspace(0.7, 1, self.PopNum)  # 生成重组概率列表，为不同的种群分配不同的重组概率
        for i in range(self.PopNum):  # 遍历种群列表
            pop = population[i]  # 得到当前种群对象
            if pop.Encoding == 'P':
                recOper = ea.Xovpmx(XOVR=Pcs[i])  # 生成部分匹配交叉算子对象
                mutOper = ea.Mutinv(Pm=float(Pms[i]))  # 生成逆转变异算子对象
            else:
                recOper = ea.Xovdp(XOVR=Pcs[i])  # 生成两点交叉算子对象
                if pop.Encoding == 'BG':
                    mutOper = ea.Mutbin(Pm=float(Pms[i]))  # 生成二进制变异算子对象
                elif pop.Encoding == 'RI':
                    mutOper = ea.Mutbga(Pm=float(Pms[i]), MutShrink=0.5, Gradient=20)  # 生成breeder GA变异算子对象
                else:
                    raise RuntimeError('编码方式必须为''BG''、''RI''或''P''.')
            self.recOpers.append(recOper)
            self.mutOpers.append(mutOper)

    def unite(self, population):
        """
        合并种群，生成联合种群。
        注：返回的unitePop不携带Field和Chrom的信息，因为其Encoding=None。
        """
        # 遍历种群列表，构造联合种群
        unitePop = ea.Population(None, None, population[0].sizes, None,  # 第一个输入参数传入None，设置Encoding为None
                                 ObjV=population[0].ObjV,
                                 FitnV=population[0].FitnV,
                                 CV=population[0].CV,
                                 Phen=population[0].Phen)
        for i in range(1, self.PopNum):
            unitePop += population[i]
        return unitePop

    def calFitness(self, population):
        """
        计算种群个体适应度，population为种群列表
        该函数直接对输入参数population中的适应度信息进行修改，因此函数不用返回任何参数。
        """
        ObjV = np.vstack(list(pop.ObjV for pop in population))
        CV = np.vstack(list(pop.CV for pop in population)) if population[0].CV is not None else population[0].CV
        FitnV = ea.scaling(ObjV, CV, self.problem.maxormins)  # 统一计算适应度
        # 为各个种群分配适应度
        idx = 0
        for i in range(self.PopNum):
            population[i].FitnV = FitnV[idx: idx + population[i].sizes]
            idx += population[i].sizes

    def EnvSelection(self, population, NUM):  # 环境选择，选择个体保留到下一代
        FitnVs = list(pop.FitnV for pop in population)
        NewChrIxs = ea.mselecting('dup', FitnVs, NUM)  # 采用基于适应度排序的直接复制选择
        for i in range(self.PopNum):
            population[i] = (population[i])[NewChrIxs[i]]
        return population

    def run(self, prophetPops=None):  # prophetPops为先知种群列表（即包含先验知识的种群列表）
        # ==========================初始化配置===========================
        self.initialization()  # 初始化算法模板的一些动态参数
        population = self.population  # 密切注意本模板的population是一个存储种群类对象的列表
        NindAll = 0  # 记录所有种群个体总数
        # ===========================准备进化============================
        for i in range(self.PopNum):  # 遍历每个种群，初始化每个种群的染色体矩阵
            NindAll += population[i].sizes
            population[i].initChrom(population[i].sizes)  # 初始化种群染色体矩阵
            self.call_aimFunc(population[i])  # 计算种群的目标函数值
            # 插入先验知识（注意：这里不会对先知种群列表prophetPops的合法性进行检查）
            if prophetPops is not None:
                population[i] = (prophetPops[i] + population[i])[:population[i].sizes]  # 插入先知种群
        self.calFitness(population)  # 统一计算适应度
        unitePop = self.unite(population)  # 得到联合种群unitePop
        # ===========================开始进化============================
        while self.terminated(unitePop) == False:
            for i in range(self.PopNum):  # 遍历种群列表，分别对各个种群进行重组和变异
                pop = population[i]  # 得到当前种群
                # 选择
                offspring = pop[ea.selecting(self.selFunc, pop.FitnV, pop.sizes)]
                # 进行进化操作
                offspring.Chrom = self.recOpers[i].do(offspring.Chrom)  # 重组
                offspring.Chrom = self.mutOpers[i].do(offspring.Encoding, offspring.Chrom, offspring.Field)  # 变异
                self.call_aimFunc(offspring)  # 计算目标函数值
                population[i] = population[i] + offspring  # 父子合并
            self.calFitness(population)  # 统一计算适应度
            population = self.EnvSelection(population, NUM=NindAll)  # 选择个体得到新一代种群
            if self.currentGen % self.migFr == 0:
                population = self.migOpers.do(population)  # 进行种群迁移
            unitePop = self.unite(population)  # 更新联合种群
        return self.finishing(unitePop)  # 调用finishing完成后续工作并返回结果

# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea  # 导入geatpy库
from sys import path as paths
from os import path

paths.append(path.split(path.split(path.realpath(__file__))[0])[0])


class moea_NSGA2_archive_templet(ea.MoeaAlgorithm):
    """
moea_NSGA2_archive_templet : class - 带全局存档的多目标进化NSGA-II算法模板
    
算法描述:
    采用带全局存档(globalNDSet)的NSGA-II进行多目标优化。
    """

    def __init__(self, problem, population):
        ea.MoeaAlgorithm.__init__(self, problem, population)  # 先调用父类构造方法
        if population.ChromNum != 1:
            raise RuntimeError('传入的种群对象必须是单染色体的种群类型。')
        self.name = 'NSGA2-archive'
        if self.problem.M < 10:
            self.ndSort = ea.ndsortESS  # 采用ENS_SS进行非支配排序
        else:
            self.ndSort = ea.ndsortTNS  # 高维目标采用T_ENS进行非支配排序，速度一般会比ENS_SS要快
        self.selFunc = 'tour'  # 选择方式，采用锦标赛选择
        if population.Encoding == 'P':
            self.recOper = ea.Xovpmx(XOVR=1)  # 生成部分匹配交叉算子对象
            self.mutOper = ea.Mutinv(Pm=1)  # 生成逆转变异算子对象
        elif population.Encoding == 'BG':
            self.recOper = ea.Xovud(XOVR=1)  # 生成均匀交叉算子对象
            self.mutOper = ea.Mutbin(Pm=None)  # 生成二进制变异算子对象，Pm设置为None时，具体数值取变异算子中Pm的默认值
        elif population.Encoding == 'RI':
            self.recOper = ea.Recsbx(XOVR=1, n=20)  # 生成模拟二进制交叉算子对象
            self.mutOper = ea.Mutpolyn(Pm=1 / self.problem.Dim, DisI=20)  # 生成多项式变异算子对象
        else:
            raise RuntimeError('编码方式必须为''BG''、''RI''或''P''.')
        self.MAXSIZE = 10 * population.sizes  # 全局非支配解存档的大小限制，默认为10倍的种群个体数

    def reinsertion(self, population, offspring, NUM, globalNDSet):

        """
        描述:
            重插入个体产生新一代种群（采用父子合并选择的策略）。
            NUM为所需要保留到下一代的个体数目，globalNDSet为全局非支配解存档。
        """

        # 父子两代合并
        population = population + offspring
        globalNDSet = population + globalNDSet  # 将population与全局存档合并
        # 非支配排序分层
        [levels, criLevel] = self.ndSort(globalNDSet.ObjV, None, None, globalNDSet.CV, self.problem.maxormins)
        # 更新全局存档
        globalNDSet = globalNDSet[np.where(levels == 1)[0]]
        if globalNDSet.CV is not None:  # CV不为None说明有设置约束条件
            globalNDSet = globalNDSet[np.where(np.all(globalNDSet.CV <= 0, 1))[0]]  # 排除非可行解
        if globalNDSet.sizes > self.MAXSIZE:
            dis = ea.crowdis(globalNDSet.ObjV, np.ones(globalNDSet.sizes))  # 计算拥挤距离
            globalNDSet = globalNDSet[np.argsort(-dis)[:self.MAXSIZE]]  # 根据拥挤距离选择符合个数限制的解保留在存档中
        # 选择个体保留到下一代
        levels = levels[: population.sizes]  # 得到与population个体对应的levels
        dis = ea.crowdis(population.ObjV, levels)  # 计算拥挤距离
        population.FitnV[:, 0] = np.argsort(np.lexsort(np.array([dis, -levels])), kind='mergesort')  # 计算适应度
        chooseFlag = ea.selecting('dup', population.FitnV, NUM)  # 调用低级选择算子dup进行基于适应度排序的选择，保留NUM个个体
        return population[chooseFlag], globalNDSet

    def run(self, prophetPop=None):  # prophetPop为先知种群（即包含先验知识的种群）
        # ==========================初始化配置===========================
        population = self.population
        NIND = population.sizes
        self.initialization()  # 初始化算法模板的一些动态参数
        # ===========================准备进化============================
        population.initChrom()  # 初始化种群染色体矩阵
        self.call_aimFunc(population)  # 计算种群的目标函数值
        # 插入先验知识（注意：这里不会对先知种群prophetPop的合法性进行检查，故应确保prophetPop是一个种群类且拥有合法的Chrom、ObjV、Phen等属性）
        if prophetPop is not None:
            population = (prophetPop + population)[:NIND]  # 插入先知种群
        [levels, criLevel] = self.ndSort(population.ObjV, NIND, None, population.CV,
                                         self.problem.maxormins)  # 对NIND个个体进行非支配分层
        population.FitnV = (1 / levels).reshape(-1, 1)  # 直接根据levels来计算初代个体的适应度
        globalNDSet = population[np.where(levels == 1)[0]]  # 创建全局存档，该全局存档贯穿进化始终，随着进化不断更新
        if globalNDSet.CV is not None:  # CV不为None说明有设置约束条件
            globalNDSet = globalNDSet[np.where(np.all(globalNDSet.CV <= 0, 1))[0]]  # 排除非可行解
        # ===========================开始进化============================
        while self.terminated(population) == False:
            # 选择个体参与进化
            offspring = population[ea.selecting(self.selFunc, population.FitnV, NIND)]
            # 对选出的个体进行进化操作
            offspring.Chrom = self.recOper.do(offspring.Chrom)  # 重组
            offspring.Chrom = self.mutOper.do(offspring.Encoding, offspring.Chrom, offspring.Field)  # 变异
            self.call_aimFunc(offspring)  # 求进化后个体的目标函数值
            population, globalNDSet = self.reinsertion(population, offspring, NIND, globalNDSet)  # 重插入生成新一代种群，同时更新全局存档
        return self.finishing(population, globalNDSet)  # 调用finishing完成后续工作并返回结果

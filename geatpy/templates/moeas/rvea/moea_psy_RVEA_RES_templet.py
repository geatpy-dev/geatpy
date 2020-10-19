# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea  # 导入geatpy库
from scipy.spatial.distance import cdist
from sys import path as paths
from os import path as path

paths.append(path.split(path.split(path.realpath(__file__))[0])[0])


class moea_psy_RVEA_RES_templet(ea.MoeaAlgorithm):
    """
moea_psy_RVEA_RES_templet : class - 带参考点再生策略的多染色体多目标进化优化RVEA算法模板(RVEA With the Reference Vector Regeneration Strategy)
    
描述:
    采用带参考点再生策略的RVEA进行多目标优化，即参考文献[1]中的RVEA*算法。
    该算法与RVEA算法类似，不过可以更好地解决具有复杂帕累托前沿面的多目标优化问题。
    该模板是内置算法模板moea_RVEA_RES_templet的多染色体版本。
    因此里面的种群对象为支持混合编码的多染色体种群类PsyPopulation类的对象。

参考文献:
    [1] Cheng R , Jin Y , Olhofer M , et al. A Reference Vector Guided 
    Evolutionary Algorithm for Many-Objective Optimization[J]. IEEE 
    Transactions on Evolutionary Computation, 2016:1-1.
    
    """

    def __init__(self, problem, population):
        ea.MoeaAlgorithm.__init__(self, problem, population)  # 先调用父类构造方法
        if population.ChromNum == 1:
            raise RuntimeError('传入的种群对象必须是多染色体的种群类型。')
        self.name = 'psy-RVEA-RES'
        self.ndSort = ea.ndsortESS  # 设置非支配排序算子
        self.selFunc = 'urs'  # 选择方式，采用无约束随机选择
        # 由于有多个染色体，因此需要用多个重组和变异算子
        self.recOpers = []
        self.mutOpers = []
        for i in range(population.ChromNum):
            if population.Encodings[i] == 'P':
                recOper = ea.Xovpmx(XOVR=1)  # 生成部分匹配交叉算子对象
                mutOper = ea.Mutinv(Pm=1)  # 生成逆转变异算子对象
            elif population.Encodings[i] == 'BG':
                recOper = ea.Xovud(XOVR=1)  # 生成均匀交叉算子对象
                mutOper = ea.Mutbin(Pm=None)  # 生成二进制变异算子对象，Pm设置为None时，具体数值取变异算子中Pm的默认值
            elif population.Encodings[i] == 'RI':
                recOper = ea.Recsbx(XOVR=1, n=20)  # 生成模拟二进制交叉算子对象
                mutOper = ea.Mutpolyn(Pm=1 / self.problem.Dim, DisI=20)  # 生成多项式变异算子对象
            else:
                raise RuntimeError('编码方式必须为''BG''、''RI''或''P''.')
            self.recOpers.append(recOper)
            self.mutOpers.append(mutOper)
        self.a = 2  # RVEA算法中的参数alpha
        self.fr = 0.1  # RVEA算法中的参数fr

    def reinsertion(self, population, offspring, refPoint):

        """
        描述:
            重插入个体产生新一代种群（采用父子合并选择的策略）。
        """

        # 父子两代合并
        population = population + offspring
        # 得到非支配个体
        [levels, criLevel] = self.ndSort(population.ObjV, None, 1, population.CV,
                                         self.problem.maxormins)  # 非支配排序，1表示只排序到第一层即非支配个体所在的层级
        population = population[np.where(levels == 1)[0]]
        # 选择个体保留到下一代
        [chooseFlag, ans] = ea.refgselect(population.ObjV, refPoint,
                                          self.problem.M * ((self.currentGen + 1) / self.MAXGEN) ** self.a,
                                          population.CV, maxormins=self.problem.maxormins)  # ans表示不使用该返回结果
        return population[chooseFlag]

    def renewRefPoint(self, ObjV, refPoint):  # 更新参考点
        _ObjV = ObjV - np.min(ObjV, 0)
        linkIdx = np.argmax(1 - cdist(_ObjV, refPoint, 'cosine'), 1)  # 找到与参考点关联的点的索引
        noLinkIdx = list(set(range(refPoint.shape[0])) - set(linkIdx))  # 找到不与参考点关联的点的索引
        refPoint[noLinkIdx, :] = np.random.rand(len(noLinkIdx), refPoint.shape[1]) * np.max(_ObjV, 0)
        return refPoint

    def run(self, prophetPop=None):  # prophetPop为先知种群（即包含先验知识的种群）
        # ==========================初始化配置===========================
        population = self.population
        self.initialization()  # 初始化算法模板的一些动态参数
        # ===========================准备进化============================
        uniformPoint, NIND = ea.crtup(self.problem.M, population.sizes)  # 生成在单位目标维度上均匀分布的参考点集
        refPoint = np.vstack([uniformPoint, np.random.rand(NIND, self.problem.M)])  # 初始化参考点（详见注释中的参考文献）
        population.initChrom(NIND)  # 初始化种群染色体矩阵，此时种群规模将调整为uniformPoint点集的大小，initChrom函数会把种群规模给重置
        self.call_aimFunc(population)  # 计算种群的目标函数值
        # 插入先验知识（注意：这里不会对先知种群prophetPop的合法性进行检查，故应确保prophetPop是一个种群类且拥有合法的Chrom、ObjV、Phen等属性）
        if prophetPop is not None:
            print('本算法需谨慎使用先验知识，有可能会导致结果比先验知识差。')
            population = (prophetPop + population)[:NIND]  # 插入先知种群
        # ===========================开始进化============================
        while self.terminated(population) == False:
            # 选择个体参与进化
            offspring = population[ea.selecting(self.selFunc, population.sizes, NIND)]
            # 进行进化操作，分别对各个种群染色体矩阵进行重组和变异
            for i in range(population.ChromNum):
                offspring.Chroms[i] = self.recOpers[i].do(offspring.Chroms[i])  # 重组
                offspring.Chroms[i] = self.mutOpers[i].do(offspring.Encodings[i], offspring.Chroms[i],
                                                          offspring.Fields[i])  # 变异
            self.call_aimFunc(offspring)  # 求进化后个体的目标函数值
            population = self.reinsertion(population, offspring, refPoint)  # 重插入生成新一代种群
            # 修改refPoint
            refPoint[NIND:, :] = self.renewRefPoint(population.ObjV, refPoint[NIND:, :])
            if (self.currentGen) % np.ceil(self.fr * self.MAXGEN) == 0:
                refPoint[:NIND, :] = uniformPoint * (np.max(population.ObjV, 0) - np.min(population.ObjV, 0))
        # 后续处理，限制种群规模（因为此时种群规模有可能大于NIND）
        [levels, criLevel] = self.ndSort(population.ObjV, NIND, None, population.CV,
                                         self.problem.maxormins)  # 对NIND个个体进行非支配分层
        population = population[
            ea.refselect(population.ObjV, levels, criLevel, NIND, uniformPoint, self.problem.maxormins)]  # 根据参考点选择个体
        return self.finishing(population)  # 调用finishing完成后续工作并返回结果

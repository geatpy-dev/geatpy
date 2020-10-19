# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea  # 导入geatpy库
from sys import path as paths
from os import path as path

paths.append(path.split(path.split(path.realpath(__file__))[0])[0])


class moea_RVEA_templet(ea.MoeaAlgorithm):
    """
moea_RVEA_templet : class - 多目标进化优化RVEA算法模板
    
算法描述:
    采用RVEA进行多目标优化。

参考文献:
    [1] Cheng R , Jin Y , Olhofer M , et al. A Reference Vector Guided 
    Evolutionary Algorithm for Many-Objective Optimization[J]. IEEE 
    Transactions on Evolutionary Computation, 2016:1-1.
    
    """

    def __init__(self, problem, population):
        ea.MoeaAlgorithm.__init__(self, problem, population)  # 先调用父类构造方法
        if population.ChromNum != 1:
            raise RuntimeError('传入的种群对象必须是单染色体的种群类型。')
        self.name = 'RVEA'
        self.selFunc = 'urs'  # 选择方式，采用无约束随机选择
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
        self.a = 2  # RVEA算法中的参数alpha
        self.fr = 0.1  # RVEA算法中的参数fr
        self.Gamma = None  # RVEA算法中的Gamma（详见参考文献的公式10），在每次更新参考点后Gamma要重置为None以便重新计算

    def reinsertion(self, population, offspring, refPoint):

        """
        描述:
            重插入个体产生新一代种群（采用父子合并选择的策略）。
        """

        # 父子两代合并
        population = population + offspring
        # 选择个体保留到下一代
        chooseFlag, self.Gamma = ea.refgselect(population.ObjV, refPoint,
                                               self.problem.M * ((self.currentGen + 1) / self.MAXGEN) ** self.a,
                                               population.CV, self.Gamma, self.problem.maxormins)
        return population[chooseFlag]

    def run(self, prophetPop=None):  # prophetPop为先知种群（即包含先验知识的种群）
        # ==========================初始化配置===========================
        population = self.population
        self.initialization()  # 初始化算法模板的一些动态参数
        # ===========================准备进化============================
        uniformPoint, NIND = ea.crtup(self.problem.M, population.sizes)  # 生成在单位目标维度上均匀分布的参考点集
        refPoint = uniformPoint.copy()  # 初始化参考点为uniformPoint
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
            # 对选出的个体进行进化操作
            offspring.Chrom = self.recOper.do(offspring.Chrom)  # 重组
            offspring.Chrom = self.mutOper.do(offspring.Encoding, offspring.Chrom, offspring.Field)  # 变异
            self.call_aimFunc(offspring)  # 求进化后个体的目标函数值
            population = self.reinsertion(population, offspring, refPoint)  # 重插入生成新一代种群
            # 修改refPoint
            if (self.currentGen) % np.ceil(self.fr * self.MAXGEN) == 0:
                refPoint = uniformPoint * (np.max(population.ObjV, 0) - np.min(population.ObjV, 0))
                self.Gamma = None  # 重置Gamma为None
        return self.finishing(population)  # 调用finishing完成后续工作并返回结果

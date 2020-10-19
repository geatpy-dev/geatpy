# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea  # 导入geatpy库
from sys import path as paths
from os import path as currentPath

paths.append(currentPath.split(currentPath.realpath(__file__))[0])
from updateNDSet import updateNDSet  # 导入该算法模板所需的外部函数


class moea_awGA_templet(ea.MoeaAlgorithm):
    """
moea_awGA_templet : class - 多目标进化优化awGA算法模板
    
算法描述:
    采用awGA进行多目标优化。
    
参考文献:
    [1] Gen M,CHeng R. Genetic Algorithms and Engineering Optimization[M]. 
    New York: John Wiley & Sons,2000
        
    """

    def __init__(self, problem, population):
        ea.MoeaAlgorithm.__init__(self, problem, population)  # 先调用父类构造方法
        if population.ChromNum != 1:
            raise RuntimeError('传入的种群对象必须是单染色体的种群类型。')
        self.name = 'awGA'
        self.selFunc = 'tour'  # 选择方式，采用锦标赛选择
        if population.Encoding == 'P':
            self.recOper = ea.Xovpmx(XOVR=1)  # 生成部分匹配交叉算子对象
            self.mutOper = ea.Mutinv(Pm=1)  # 生成逆转变异算子对象
        elif population.Encoding == 'BG':
            self.recOper = ea.Xovud(XOVR=1)  # 生成部均匀交叉算子对象
            self.mutOper = ea.Mutbin(Pm=None)  # 生成二进制变异算子对象，Pm设置为None时，具体数值取变异算子中Pm的默认值
        elif population.Encoding == 'RI':
            self.recOper = ea.Xovud(XOVR=1)  # 生成部均匀交叉算子对象
            self.mutOper = ea.Mutuni(Pm=1 / self.problem.Dim, Alpha=False, Middle=False)  # 生成均匀变异算子对象
            self.extraMutOper = ea.Mutgau(Pm=1 / self.problem.Dim, Sigma3=False,
                                          Middle=False)  # 额外生成一个高斯变异算子对象，对标准差放大3倍
        else:
            raise RuntimeError('编码方式必须为''BG''、''RI''或''P''.')
        self.MAXSIZE = population.sizes  # 非支配解集大小限制

    def run(self, prophetPop=None):  # prophetPop为先知种群（即包含先验知识的种群）
        # ==========================初始化配置===========================
        problem = self.problem
        population = self.population
        NIND = population.sizes
        MAXSIZE = self.MAXSIZE
        if MAXSIZE is None:  # 检查MAXSIZE，默认取2倍的种群规模
            MAXSIZE = 2 * NIND
        self.initialization()  # 初始化算法模板的一些动态参数
        # ===========================准备进化============================
        population.initChrom(NIND)  # 初始化种群染色体矩阵
        self.call_aimFunc(population)  # 计算种群的目标函数值
        NDSet = updateNDSet(population, problem.maxormins, MAXSIZE)  # 计算适应度和得到全局非支配种群
        self.evalsNum = population.sizes  # 记录评价次数
        # 插入先验知识（注意：这里不会对先知种群prophetPop的合法性进行检查，故应确保prophetPop是一个种群类且拥有合法的Chrom、ObjV、Phen等属性）
        if prophetPop is not None:
            population = (prophetPop + population)[:NIND]  # 插入先知种群
        # ===========================开始进化============================
        while self.terminated(population) == False:
            uniChrom = np.unique(NDSet.Chrom, axis=0)
            repRate = 1 - uniChrom.shape[0] / NDSet.sizes  # 计算NDSet中的重复率
            # 选择个体去进化形成子代
            offspring = population[ea.selecting(self.selFunc, population.FitnV, NIND)]
            offspring.Chrom = self.recOper.do(offspring.Chrom)  # 重组
            offspring.Chrom = self.mutOper.do(offspring.Encoding, offspring.Chrom, offspring.Field)  # 变异
            if population.Encoding == 'RI' and repRate > 0.1:
                offspring.Chrom = self.extraMutOper.do(offspring.Encoding, offspring.Chrom, offspring.Field)  # 执行额外的变异
            self.call_aimFunc(offspring)  # 求进化后个体的目标函数值
            population = population + offspring  # 父代种群和育种种群合并
            NDSet = updateNDSet(population, problem.maxormins, MAXSIZE, NDSet)  # 计算合并种群的适应度及更新NDSet
            # 保留个体到下一代
            population = population[ea.selecting('dup', population.FitnV, NIND)]  # 选择，保留NIND个个体
        if NDSet.CV is not None:  # CV不为None说明有设置约束条件
            NDSet = NDSet[np.where(np.all(NDSet.CV <= 0, 1))[0]]  # 最后要彻底排除非可行解
        return self.finishing(population, NDSet)  # 调用finishing完成后续工作并返回结果

# -*- coding: utf-8 -*-
import geatpy as ea  # 导入geatpy库
from sys import path as paths
from os import path as path

paths.append(path.split(path.split(path.realpath(__file__))[0])[0])


class moea_NSGA3_templet(ea.MoeaAlgorithm):
    """
moea_NSGA3_templet : class - 多目标进化优化NSGA-III算法模板
    
算法描述:
    采用NSGA-III进行多目标优化。
    注意：在初始化染色体时，种群规模会被修正为NSGA-III所用的参考点集的大小。

参考文献:
    [1] Deb K , Jain H . An Evolutionary Many-Objective Optimization Algorithm 
    Using Reference-Point-Based Nondominated Sorting Approach, Part I: 
    Solving Problems With Box Constraints[J]. IEEE Transactions on 
    Evolutionary Computation, 2014, 18(4):577-601.
    
    """

    def __init__(self, problem, population):
        ea.MoeaAlgorithm.__init__(self, problem, population)  # 先调用父类构造方法
        if population.ChromNum != 1:
            raise RuntimeError('传入的种群对象必须是单染色体的种群类型。')
        self.name = 'NSGA3'
        if self.problem.M < 10:
            self.ndSort = ea.ndsortESS  # 采用ENS_SS进行非支配排序
        else:
            self.ndSort = ea.ndsortTNS  # 高维目标采用T_ENS进行非支配排序，速度一般会比ENS_SS要快
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

    def reinsertion(self, population, offspring, NUM, uniformPoint):

        """
        描述:
            重插入个体产生新一代种群（采用父子合并选择的策略）。
            NUM为所需要保留到下一代的个体数目。
        """

        # 父子两代合并
        population = population + offspring
        # 选择个体保留到下一代
        [levels, criLevel] = self.ndSort(population.ObjV, NUM, None, population.CV,
                                         self.problem.maxormins)  # 对NUM个个体进行非支配分层
        chooseFlag = ea.refselect(population.ObjV, levels, criLevel, NUM, uniformPoint,
                                  self.problem.maxormins)  # 根据参考点的“入龛”个体筛选
        return population[chooseFlag]

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
        # ===========================开始进化============================
        while self.terminated(population) == False:
            # 选择个体参与进化
            offspring = population[ea.selecting(self.selFunc, population.sizes, NIND)]
            # 对选出的个体进行进化操作
            offspring.Chrom = self.recOper.do(offspring.Chrom)  # 重组
            offspring.Chrom = self.mutOper.do(offspring.Encoding, offspring.Chrom, offspring.Field)  # 变异
            self.call_aimFunc(offspring)  # 求进化后个体的目标函数值
            # 重插入生成新一代种群
            population = self.reinsertion(population, offspring, NIND, uniformPoint)
        return self.finishing(population)  # 调用finishing完成后续工作并返回结果

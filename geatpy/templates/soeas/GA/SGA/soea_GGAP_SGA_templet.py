# -*- coding: utf-8 -*-
import geatpy as ea  # 导入geatpy库
import numpy as np
from sys import path as paths
from os import path

paths.append(path.split(path.split(path.realpath(__file__))[0])[0])


class soea_GGAP_SGA_templet(ea.SoeaAlgorithm):
    """
soea_GGAP_SGA_templet : class - Generational Gap Simple GA templet(带代沟的简单遗传算法模板)

算法描述:
    本模板实现的是带代沟的简单遗传算法，
    它在SGA算法模板的基础上增加“代沟”，用于控制使用多少个子代替换父代来形成新一代种群，算法流程如下：
    1) 根据编码规则初始化N个个体的种群。
    2) 若满足停止条件则停止，否则继续执行。
    3) 对当前种群进行统计分析，比如记录其最优个体、平均适应度等等。
    4) 独立地从当前种群中选取N个母体。
    5) 独立地对这N个母体进行交叉操作。
    6) 独立地对这N个交叉后的个体进行变异，并根据代沟从中选择N'个个体替换父代最差的N'个个体，得到下一代种群。
    7) 回到第2步。
    
"""

    def __init__(self, problem, population):
        ea.SoeaAlgorithm.__init__(self, problem, population)  # 先调用父类构造方法
        if population.ChromNum != 1:
            raise RuntimeError('传入的种群对象必须是单染色体的种群类型。')
        self.name = 'GGAP-SGA'
        self.selFunc = 'rws'  # 轮盘赌选择算子
        if population.Encoding == 'P':
            self.recOper = ea.Xovpmx(XOVR=0.7)  # 生成部分匹配交叉算子对象
            self.mutOper = ea.Mutinv(Pm=0.5)  # 生成逆转变异算子对象
        else:
            self.recOper = ea.Xovdp(XOVR=0.7)  # 生成两点交叉算子对象
            if population.Encoding == 'BG':
                self.mutOper = ea.Mutbin(Pm=None)  # 生成二进制变异算子对象，Pm设置为None时，具体数值取变异算子中Pm的默认值
            elif population.Encoding == 'RI':
                self.mutOper = ea.Mutbga(Pm=1 / self.problem.Dim, MutShrink=0.5, Gradient=20)  # 生成breeder GA变异算子对象
            else:
                raise RuntimeError('编码方式必须为''BG''、''RI''或''P''.')
        self.GGAP = 0.9  # 代沟，表示使用多少比例的子代替换父代来形成新一代种群

    def reinsertion(self, population, offspring, GGAP_NUM):

        """ 重插入 """
        replaceIdx = np.argsort(population.FitnV.T[0])[:GGAP_NUM].astype(int)  # 计算父代中要被替换的个体索引
        insertIdx = np.argsort(-offspring.FitnV.T[0])[:GGAP_NUM].astype(int)  # 计算子代中需要选择进行重插入的个体索引
        population[replaceIdx] = offspring[insertIdx]
        return population

    def run(self, prophetPop=None):  # prophetPop为先知种群（即包含先验知识的种群）
        # ==========================初始化配置===========================
        population = self.population
        NIND = population.sizes
        GGAP_NUM = int(np.ceil(NIND * self.GGAP))  # 计算每一代替换个体的个数
        self.initialization()  # 初始化算法模板的一些动态参数
        # ===========================准备进化============================
        population.initChrom(NIND)  # 初始化种群染色体矩阵
        self.call_aimFunc(population)  # 计算种群的目标函数值
        # 插入先验知识（注意：这里不会对先知种群prophetPop的合法性进行检查，故应确保prophetPop是一个种群类且拥有合法的Chrom、ObjV、Phen等属性）
        if prophetPop is not None:
            population = (prophetPop + population)[:NIND]  # 插入先知种群
        population.FitnV = ea.scaling(population.ObjV, population.CV, self.problem.maxormins)  # 计算适应度
        # ===========================开始进化============================
        while self.terminated(population) == False:
            # 选择
            offspring = population[ea.selecting(self.selFunc, population.FitnV, NIND)]
            # 进行进化操作
            offspring.Chrom = self.recOper.do(offspring.Chrom)  # 重组
            offspring.Chrom = self.mutOper.do(offspring.Encoding, offspring.Chrom, offspring.Field)  # 变异
            self.call_aimFunc(offspring)  # 计算目标函数值
            offspring.FitnV = ea.scaling(offspring.ObjV, offspring.CV, self.problem.maxormins)  # 计算适应度
            # 根据代沟把子代重插入到父代生成新一代种群
            population = self.reinsertion(population, offspring, GGAP_NUM)
            population.FitnV = ea.scaling(population.ObjV, population.CV, self.problem.maxormins)  # 计算适应度
        return self.finishing(population)  # 调用finishing完成后续工作并返回结果

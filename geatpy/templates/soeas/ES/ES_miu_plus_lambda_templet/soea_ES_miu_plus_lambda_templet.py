# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea  # 导入geatpy库
from sys import path as paths
from os import path

paths.append(path.split(path.split(path.realpath(__file__))[0])[0])


class soea_ES_miu_plus_lambda_templet(ea.SoeaAlgorithm):
    """
soea_ES_miu_plus_lambda_templet : class - (μ+λ)进化策略模板

算法描述:
    本模板实现的是(μ+λ)进化策略[1]。

参考文献:
    [1] Beyer H G , Schwefel H P . Evolution strategies – A comprehensive 
    introduction[J]. Natural Computing, 2002, 1(1):3-52.

"""

    def __init__(self, problem, population):
        ea.SoeaAlgorithm.__init__(self, problem, population)  # 先调用父类构造方法
        if population.ChromNum != 1:
            raise RuntimeError('传入的种群对象必须是单染色体的种群类型。')
        self.name = '(μ+λ)ES'
        if population.Encoding != 'RI':
            raise RuntimeError('编码方式必须为''RI''.')
        self.selFunc = 'urs'  # 无约束的随机选择算子
        self.NSel = int(0.5 * population.sizes)  # 这里用NSel代指算法中的lambda，默认设为种群规模的0.5倍
        self.recOper = ea.Recint(RecOpt=None, Half_N=self.NSel, Alpha=0.5)  # 默认采用全局重组方式的中间重组

    def run(self, prophetPop=None):  # prophetPop为先知种群（即包含先验知识的种群）
        # ==========================初始化配置===========================
        population = self.population
        NIND = population.sizes
        self.initialization()  # 初始化算法模板的一些动态参数
        # ===========================准备进化============================
        population.initChrom(NIND)  # 初始化种群染色体矩阵
        self.call_aimFunc(population)  # 计算种群的目标函数值
        # 插入先验知识（注意：这里不会对先知种群prophetPop的合法性进行检查，故应确保prophetPop是一个种群类且拥有合法的Chrom、ObjV、Phen等属性）
        if prophetPop is not None:
            population = (prophetPop + population)[:NIND]  # 插入先知种群
        population.FitnV = ea.scaling(population.ObjV, population.CV, self.problem.maxormins)  # 计算适应度
        Sigma3 = np.random.rand(population.sizes, population.Lind) * (
                    population.Field[1, :] - population.Field[0, :]) * 0.5  # 初始化高斯变异算子的Sigma3
        # ===========================开始进化============================
        while self.terminated(population) == False:
            # 选择
            choose_index = ea.selecting(self.selFunc, population.FitnV, self.NSel)
            offspring = population[choose_index]
            offspring_Sigma3 = Sigma3[choose_index]
            # 进行进化操作
            offspring.Chrom = self.recOper.do(offspring.Chrom)  # 重组
            offspring_Sigma3 = self.recOper.do(offspring_Sigma3)  # 对offspring_Sigma3进行重组
            for i in range(len(choose_index)):
                offspring.Chrom[i, :] = ea.mutgau(offspring.Encoding, offspring.Chrom[i, :].reshape(1, -1),
                                                  offspring.Field, 1, offspring_Sigma3[i]).reshape(-1)  # 高斯变异
            self.call_aimFunc(offspring)  # 计算目标函数值
            population = population + offspring  # 父子合并
            Sigma3 = np.vstack([Sigma3, offspring_Sigma3])
            population.FitnV = ea.scaling(population.ObjV, population.CV, self.problem.maxormins)  # 计算适应度
            # 得到新一代种群
            choose_index = ea.selecting('dup', population.FitnV, NIND)  # 采用基于适应度排序的直接复制选择
            population = population[choose_index]
            Sigma3 = Sigma3[choose_index]
        return self.finishing(population)  # 调用finishing完成后续工作并返回结果

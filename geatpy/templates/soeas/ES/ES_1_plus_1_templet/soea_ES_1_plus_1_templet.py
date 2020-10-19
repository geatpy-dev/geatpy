# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea  # 导入geatpy库
from sys import path as paths
from os import path

paths.append(path.split(path.split(path.realpath(__file__))[0])[0])


class soea_ES_1_plus_1_templet(ea.SoeaAlgorithm):
    """
soea_ES_1_plus_1_templet : class - (1+1)进化策略模板

算法描述:
    本模板实现的是(1+1)进化策略。算法流程如下：
    1) 根据编码规则初始化N个个体的种群。
    2) 若满足停止条件则停止，否则继续执行。
    3) 对当前种群进行统计分析，比如记录其最优个体、平均适应度等等。
    4) 初始化控制高斯变异中的标准差Sigma(Geatpy中的高斯变异算子传入的是3倍的标准差即Sigma3)。
    5) 独立地对这种群个体进行高斯变异，得到试验种群。
    6) 在当前种群和实验种群之间采用一对一生存者选择方法得到新一代种群，
    同时统计新一代种群中有多少个个体继承自实验种群（即变异成功率）。
    7) 根据变异成功率修改Sigma3。
    8) 回到第2步。

参考文献:
    [1] Beyer H G , Schwefel H P . Evolution strategies – A comprehensive 
    introduction[J]. Natural Computing, 2002, 1(1):3-52.

"""

    def __init__(self, problem, population):
        ea.SoeaAlgorithm.__init__(self, problem, population)  # 先调用父类构造方法
        if population.ChromNum != 1:
            raise RuntimeError('传入的种群对象必须是单染色体的种群类型。')
        self.name = '(1+1)ES'
        if population.Encoding != 'RI':
            raise RuntimeError('编码方式必须为''RI''.')

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
        Sigma3 = 0.5 * (population.Field[1, :] - population.Field[0, :])  # 初始化高斯变异算子的Sigma3
        # ===========================开始进化============================
        while self.terminated(population) == False:
            # 进行进化操作
            experimentPop = ea.Population(population.Encoding, population.Field, NIND)  # 存储试验种群
            experimentPop.Chrom = ea.mutgau(population.Encoding, population.Chrom, population.Field, 1, Sigma3)  # 高斯变异
            self.call_aimFunc(experimentPop)  # 计算目标函数值
            tempPop = population + experimentPop  # 临时合并，以调用otos进行一对一生存者选择
            tempPop.FitnV = ea.scaling(tempPop.ObjV, tempPop.CV, self.problem.maxormins)  # 计算适应度
            chooseIdx = ea.selecting('otos', tempPop.FitnV, NIND)  # 采用One-to-One Survivor选择
            population = tempPop[chooseIdx]  # 产生新一代种群
            # 利用1/5规则调整变异压缩概率
            successfulRate = len(np.where(chooseIdx >= NIND)[0]) / (2 * NIND)
            if successfulRate < 1 / 5:
                Sigma3 *= 0.817
            elif successfulRate > 1 / 5:
                Sigma3 /= 0.817
        return self.finishing(population)  # 调用finishing完成后续工作并返回结果

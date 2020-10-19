# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea  # 导入geatpy库
from sys import path as paths
from os import path

paths.append(path.split(path.split(path.realpath(__file__))[0])[0])


class soea_studGA_templet(ea.SoeaAlgorithm):
    """
soea_studGA_templet.py - Stud GA templet(种马遗传算法模板)

算法描述:
    本模板实现的是种马遗传算法。算法流程详见参考文献[1]。

参考文献:
    [1] Khatib W , Fleming P J . The stud GA: A mini revolution?[C]// International 
    Conference on Parallel Problem Solving from Nature. Springer, Berlin, Heidelberg, 1998.
    
"""

    def __init__(self, problem, population):
        ea.SoeaAlgorithm.__init__(self, problem, population)  # 先调用父类构造方法
        if population.ChromNum != 1:
            raise RuntimeError('传入的种群对象必须是单染色体的种群类型。')
        self.name = 'studGA'
        self.selFunc = 'tour'  # 锦标赛选择算子
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
        # ===========================开始进化============================
        while self.terminated(population) == False:
            bestIdx = np.argmax(population.FitnV, axis=0)  # 得到当代的最优个体的索引, 设置axis=0可使得返回一个向量
            studPop = population[np.tile(bestIdx, (NIND // 2))]  # 复制最优个体NIND//2份，组成一个“种马种群”
            restPop = population[np.where(np.arange(NIND) != bestIdx)[0]]  # 得到除去精英个体外其它个体组成的种群
            # 选择个体，以便后面与种马种群进行交配
            tempPop = restPop[ea.selecting(self.selFunc, restPop.FitnV, (NIND - studPop.sizes))]
            # 将种马种群与选择出来的个体进行合并
            population = studPop + tempPop
            # 进行进化操作
            population.Chrom = self.recOper.do(population.Chrom)  # 重组
            population.Chrom = self.mutOper.do(population.Encoding, population.Chrom, population.Field)  # 变异
            self.call_aimFunc(population)
            population.FitnV = ea.scaling(population.ObjV, population.CV, self.problem.maxormins)  # 计算适应度
        return self.finishing(population)  # 调用finishing完成后续工作并返回结果

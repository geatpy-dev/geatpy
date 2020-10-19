# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea  # 导入geatpy库
from sys import path as paths
from os import path

paths.append(path.split(path.split(path.realpath(__file__))[0])[0])


class soea_psy_studGA_templet(ea.SoeaAlgorithm):
    """
soea_psy_studGA_templet.py - Polysomy Stud GA templet(多染色体种马遗传算法模板)

模板说明:
    该模板是内置算法模板soea_studGA_templet的多染色体版本，
    因此里面的种群对象为支持混合编码的多染色体种群类PsyPopulation类的对象。
    
算法描述:
    本模板实现的是种马遗传算法。算法流程详见参考文献[1]。

模板使用注意:
    本模板调用的目标函数形如：aimFunc(pop), 
    其中pop为种群类的对象，代表一个种群，
    pop对象的Phen属性（即种群染色体的表现型）等价于种群所有个体的决策变量组成的矩阵，
    该函数根据该Phen计算得到种群所有个体的目标函数值组成的矩阵，并将其赋值给pop对象的ObjV属性。
    若有约束条件，则在计算违反约束程度矩阵CV后赋值给pop对象的CV属性（详见Geatpy数据结构）。
    该函数不返回任何的返回值，求得的目标函数值保存在种群对象的ObjV属性中，
                          违反约束程度矩阵保存在种群对象的CV属性中。
    例如：population为一个种群对象，则调用aimFunc(population)即可完成目标函数值的计算，
         此时可通过population.ObjV得到求得的目标函数值，population.CV得到违反约束程度矩阵。
    若不符合上述规范，则请修改算法模板或自定义新算法模板。

参考文献:
    [1] Khatib W , Fleming P J . The stud GA: A mini revolution?[C]// International 
    Conference on Parallel Problem Solving from Nature. Springer, Berlin, Heidelberg, 1998.
    
"""

    def __init__(self, problem, population):
        ea.SoeaAlgorithm.__init__(self, problem, population)  # 先调用父类构造方法
        if population.ChromNum == 1:
            raise RuntimeError('传入的种群对象必须是多染色体的种群类型。')
        self.name = 'psy-studGA'
        self.problem = problem
        self.population = population
        self.selFunc = 'tour'  # 锦标赛选择算子
        # 由于有多个染色体，因此需要用多个重组和变异算子
        self.recOpers = []
        self.mutOpers = []
        for i in range(population.ChromNum):
            if population.Encodings[i] == 'P':
                recOper = ea.Xovpmx(XOVR=0.7)  # 生成部分匹配交叉算子对象
                mutOper = ea.Mutinv(Pm=0.5)  # 生成逆转变异算子对象
            else:
                recOper = ea.Xovdp(XOVR=0.7)  # 生成两点交叉算子对象
                if population.Encodings[i] == 'BG':
                    mutOper = ea.Mutbin(Pm=None)  # 生成二进制变异算子对象，Pm设置为None时，具体数值取变异算子中Pm的默认值
                elif population.Encodings[i] == 'RI':
                    mutOper = ea.Mutbga(Pm=1 / self.problem.Dim, MutShrink=0.5, Gradient=20)  # 生成breeder GA变异算子对象
                else:
                    raise RuntimeError('编码方式必须为''BG''、''RI''或''P''.')
            self.recOpers.append(recOper)
            self.mutOpers.append(mutOper)

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
            # 进行进化操作，分别对各种编码的染色体进行重组和变异
            for i in range(population.ChromNum):
                population.Chroms[i] = self.recOpers[i].do(population.Chroms[i])  # 重组
                population.Chroms[i] = self.mutOpers[i].do(population.Encodings[i], population.Chroms[i],
                                                           population.Fields[i])  # 变异
            self.call_aimFunc(population)
            population.FitnV = ea.scaling(population.ObjV, population.CV, self.problem.maxormins)  # 计算适应度
        return self.finishing(population)  # 调用finishing完成后续工作并返回结果

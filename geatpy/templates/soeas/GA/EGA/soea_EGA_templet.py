# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea # 导入geatpy库
from sys import path as paths
from os import path
paths.append(path.split(path.split(path.realpath(__file__))[0])[0])

class soea_EGA_templet(ea.SoeaAlgorithm):
    
    """
soea_EGA_templet.py - Elitist Reservation GA templet(精英保留的遗传算法模板)

算法描述:
    本模板实现的是基于杰出保留的单目标遗传算法。算法流程如下：
    1) 根据编码规则初始化N个个体的种群。
    2) 若满足停止条件则停止，否则继续执行。
    3) 对当前种群进行统计分析，比如记录其最优个体、平均适应度等等。
    4) 独立地从当前种群中选取N-1个母体。
    5) 独立地对这N-1个母体进行交叉操作。
    6) 独立地对这N-1个交叉后的个体进行变异。
    7) 计算当代种群的最优个体，并把它插入到这N-1个交叉后的个体的第一位，得到新一代种群。
    8) 回到第2步。

模板使用注意:
    本模板调用的目标函数形如：aimFunc(pop), 
    其中pop为Population类的对象，代表一个种群，
    pop对象的Phen属性（即种群染色体的表现型）等价于种群所有个体的决策变量组成的矩阵，
    该函数根据该Phen计算得到种群所有个体的目标函数值组成的矩阵，并将其赋值给pop对象的ObjV属性。
    若有约束条件，则在计算违反约束程度矩阵CV后赋值给pop对象的CV属性（详见Geatpy数据结构）。
    该函数不返回任何的返回值，求得的目标函数值保存在种群对象的ObjV属性中，
                          违反约束程度矩阵保存在种群对象的CV属性中。
    例如：population为一个种群对象，则调用aimFunc(population)即可完成目标函数值的计算，
         此时可通过population.ObjV得到求得的目标函数值，population.CV得到违反约束程度矩阵。
    若不符合上述规范，则请修改算法模板或自定义新算法模板。
    
"""
    
    def __init__(self, problem, population):
        ea.SoeaAlgorithm.__init__(self, problem, population) # 先调用父类构造方法
        if str(type(population)) != "<class 'Population.Population'>":
            raise RuntimeError('传入的种群对象必须为Population类型')
        self.name = 'EGA'
        self.selFunc = 'tour' # 锦标赛选择算子
        if population.Encoding == 'P':
            self.recFunc = 'xovpmx' # 部分匹配交叉
            self.mutFunc = 'mutinv' # 染色体片段逆转变异
        else:
            self.recFunc = 'xovdp'  # 两点交叉
            if population.Encoding == 'BG':
                self.mutFunc = 'mutbin' # 二进制变异
            elif population.Encoding == 'RI':
                self.mutFunc = 'mutbga' # breeder GA中的变异算子
            else:
                raise RuntimeError('编码方式必须为''BG''、''RI''或''P''.')
        self.pc = 1 # 重组概率
        self.pm = 1 # 整条染色体的变异概率
        
    def run(self):
        #==========================初始化配置===========================
        population = self.population
        NIND = population.sizes
        self.initialization() # 初始化算法模板的一些动态参数
        #===========================准备进化============================
        if population.Chrom is None:
            population.initChrom(NIND) # 初始化种群染色体矩阵（内含染色体解码，详见Population类的源码）
        self.problem.aimFunc(population) # 计算种群的目标函数值
        population.FitnV = ea.scaling(self.problem.maxormins * population.ObjV, population.CV) # 计算适应度
        self.evalsNum = population.sizes # 记录评价次数
        #===========================开始进化============================
        while self.terminated(population) == False:
            bestIndi = population[np.argmax(population.FitnV, 0)] # 得到当代的最优个体
            # 选择
            offspring = population[ea.selecting(self.selFunc, population.FitnV, NIND - 1)]
            # 进行进化操作
            offspring.Chrom = ea.recombin(self.recFunc, offspring.Chrom, self.pc) # 重组
            offspring.Chrom = ea.mutate(self.mutFunc, offspring.Encoding, offspring.Chrom, offspring.Field, self.pm) # 变异
            # 求进化后个体的目标函数值
            offspring.Phen = offspring.decoding() # 染色体解码
            self.problem.aimFunc(offspring) # 计算目标函数值
            self.evalsNum += offspring.sizes # 更新评价次数
            population = bestIndi + offspring # 更新种群
            population.FitnV = ea.scaling(self.problem.maxormins * population.ObjV, population.CV) # 计算适应度
        
        return self.finishing(population) # 调用finishing完成后续工作并返回结果
    
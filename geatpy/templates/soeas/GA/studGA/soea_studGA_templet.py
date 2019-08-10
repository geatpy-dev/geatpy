# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea # 导入geatpy库
from sys import path as paths
from os import path
paths.append(path.split(path.split(path.realpath(__file__))[0])[0])

class soea_studGA_templet(ea.SoeaAlgorithm):
    
    """
soea_studGA_templet.py - Stud GA templet(种马遗传算法模板)

算法描述:
    本模板实现的是种马遗传算法。算法流程详见参考文献[1]。

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

参考文献:
    [1] Khatib W , Fleming P J . The stud GA: A mini revolution?[C]// International 
    Conference on Parallel Problem Solving from Nature. Springer, Berlin, Heidelberg, 1998.
    
"""
    
    def __init__(self, problem, population):
        ea.SoeaAlgorithm.__init__(self, problem, population) # 先调用父类构造方法
        if str(type(population)) != "<class 'Population.Population'>":
            raise RuntimeError('传入的种群对象必须为Population类型')
        self.name = 'studGA'
        self.problem = problem
        self.population = population
        self.selFunc = 'tour' # 锦标赛选择算子
        if population.Encoding == 'P':
            self.recFunc = 'xovpmx' # 部分匹配交叉
            self.mutFunc = 'mutinv' # 染色体片段逆转变异
        else:
            self.recFunc = 'xovdp' # 两点交叉
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
            bestIdx = np.argmax(population.FitnV, axis = 0) # 得到当代的最优个体的索引, 设置axis=0可使得返回一个向量
            studPop = population[np.tile(bestIdx, (NIND//2))] # 复制最优个体NIND//2份，组成一个“种马种群”
            restPop = population[np.where(np.array(range(NIND)) != bestIdx)[0]] # 得到除去精英个体外其它个体组成的种群
            # 选择个体，以便后面与种马种群进行交配
            tempPop = restPop[ea.selecting(self.selFunc, restPop.FitnV, (NIND - studPop.sizes))]
            # 将种马种群与选择出来的个体进行合并
            population = studPop + tempPop
            # 进行进化操作
            population.Chrom = ea.recombin(self.recFunc, population.Chrom, self.pc) # 重组
            population.Chrom = ea.mutate(self.mutFunc, population.Encoding, population.Chrom, population.Field, self.pm) # 变异
            # 求进化后个体的目标函数值
            population.Phen = population.decoding() # 染色体解码
            self.problem.aimFunc(population)
            self.evalsNum += population.sizes # 更新评价次数
            population.FitnV = ea.scaling(self.problem.maxormins * population.ObjV, population.CV) # 计算适应度
        
        return self.finishing(population) # 调用finishing完成后续工作并返回结果
    
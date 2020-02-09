# -*- coding: utf-8 -*-
import time
import numpy as np
import geatpy as ea # 导入geatpy库
from sys import path as paths
from os import path as currentPath
paths.append(currentPath.split(currentPath.realpath(__file__))[0])
from updateNDSet import updateNDSet # 导入该算法模板所需的外部函数

class moea_awGA_templet(ea.MoeaAlgorithm):
    
    """
moea_awGA_templet : class - 多目标进化优化awGA算法模板
    
算法描述:
    采用awGA进行多目标优化。

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
    [1] Gen M,CHeng R. Genetic Algorithms and Engineering Optimization[M]. 
    New York: John Wiley & Sons,2000
        
    """
    
    def __init__(self, problem, population):
        ea.MoeaAlgorithm.__init__(self, problem, population) # 先调用父类构造方法
        if str(type(population)) != "<class 'Population.Population'>":
            raise RuntimeError('传入的种群对象必须为Population类型')
        self.name = 'awGA'
        self.selFunc = 'tour' # 选择方式，采用锦标赛选择
        if population.Encoding == 'P':
            self.recOper = ea.Xovpmx(XOVR = 1) # 生成部分匹配交叉算子对象
            self.mutOper = ea.Mutinv(Pm = 1) # 生成逆转变异算子对象
        elif population.Encoding == 'BG':
            self.recOper = ea.Xovud(XOVR = 1) # 生成部均匀交叉算子对象
            self.mutOper = ea.Mutbin(Pm = None) # 生成二进制变异算子对象，Pm设置为None时，具体数值取变异算子中Pm的默认值
        elif population.Encoding == 'RI':
            self.recOper = ea.Xovud(XOVR = 1) # 生成部均匀交叉算子对象
            self.mutOper = ea.Mutuni(Pm = 1/self.problem.Dim, Alpha = False, MutShrink = 1, Middle = False) # 生成均匀变异算子对象
            self.extraMutOper = ea.Mutgau(Pm = 1/self.problem.Dim, Sigma3 = False, MutShrink = 3, Middle = False) # 额外生成一个高斯变异算子对象，对标准差放大3倍
        else:
            raise RuntimeError('编码方式必须为''BG''、''RI''或''P''.')
            
        self.MAXSIZE = population.sizes # 非支配解集大小限制

    def run(self, prophetPop = None): # prophetPop为先知种群（即包含先验知识的种群）
        #==========================初始化配置===========================
        problem = self.problem
        population = self.population
        NIND = population.sizes
        MAXSIZE = self.MAXSIZE
        if MAXSIZE is None: # 检查MAXSIZE，默认取2倍的种群规模
            MAXSIZE = 2 * NIND
        self.initialization() # 初始化算法模板的一些动态参数
        #===========================准备进化============================
        if population.Chrom is None:
            population.initChrom(NIND) # 初始化种群染色体矩阵（内含解码，详见Population类的源码）
        else:
            population.Phen = population.decoding() # 染色体解码
        self.problem.aimFunc(population) # 计算种群的目标函数值
        NDSet = updateNDSet(population, problem.maxormins, MAXSIZE) # 计算适应度和得到全局非支配种群
        self.evalsNum = population.sizes # 记录评价次数
        # 插入先验知识（注意：这里不会对先知种群prophetPop的合法性进行检查，故应确保prophetPop是一个种群类且拥有合法的Chrom、ObjV、Phen等属性）
        if prophetPop is not None:
            population = (prophetPop + population)[:NIND] # 插入先知种群
        #===========================开始进化============================
        while self.terminated(population) == False:
            uniChrom = np.unique(NDSet.Chrom, axis = 0)
            repRate = 1 - uniChrom.shape[0] / NDSet.sizes # 计算NDSet中的重复率
            # 选择个体去进化形成子代
            offspring = population[ea.selecting(self.selFunc, population.FitnV, NIND)]
            offspring.Chrom = self.recOper.do(offspring.Chrom) # 重组
            offspring.Chrom = self.mutOper.do(offspring.Encoding, offspring.Chrom, offspring.Field) # 变异
            if population.Encoding == 'RI' and repRate > 0.1:
                offspring.Chrom = self.extraMutOper.do(offspring.Encoding, offspring.Chrom, offspring.Field) # 执行额外的变异
            offspring.Phen = offspring.decoding() # 染色体解码
            self.problem.aimFunc(offspring) # 求进化后个体的目标函数值
            self.evalsNum += offspring.sizes # 更新评价次数
            # 父代种群和育种种群合并
            population = population + offspring
            NDSet = updateNDSet(population, problem.maxormins, MAXSIZE, NDSet) # 计算合并种群的适应度及更新NDSet
            # 保留个体到下一代
            population = population[ea.selecting('dup', population.FitnV, NIND)] # 选择，保留NIND个个体
        NDSet = NDSet[np.where(np.all(NDSet.CV <= 0, 1))[0]] # 最后要彻底排除非可行解
        self.passTime += time.time() - self.timeSlot # 更新用时记录
        #=========================绘图及输出结果=========================
        if self.drawing != 0:
            ea.moeaplot(NDSet.ObjV, Label = 'Pareto Front', saveFlag = True, gridFlag = True)
        # 返回帕累托最优集
        return NDSet
    
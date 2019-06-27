# -*- coding: utf-8 -*-
import time
import numpy as np
import geatpy as ea # 导入geatpy库
from sys import path as paths
from os import path as currentPath
paths.append(currentPath.split(currentPath.realpath(__file__))[0])
from updateNDSet import updateNDSet # 导入该算法模板所需的外部函数

class moea_awGA_templet(ea.Algorithm):
    
    """
moea_awGA_templet : class - 基于awGA算法求解多目标优化问题的进化算法模板
    
算法描述:
    采用awGA进行多目标优化
    
参考文献:
    [1] Gen M,CHeng R. Genetic Algorithms and Engineering Optimization[M]. 
    New York: John Wiley & Sons,2000
    
模板使用注意：
    1.本模板调用的目标函数形如：[ObjV,CV] = aimFuc(Vars, CV), 
      其中Vars表示决策变量矩阵, CV为种群的违反约束程度矩阵(详见Geatpy数据结构),
    若不符合上述规范，则请修改算法模板或自定义新算法模板
        
    """
    
    def __init__(self, problem, population):
        self.name = 'awGA'
        self.problem = problem
        self.population = population
        self.selFunc = 'tour' # 选择方式，采用锦标赛选择
        if population.Encoding == 'B' or population.Encoding == 'G':
            self.recFunc = 'xovud' # 均匀交叉
            self.mutFunc = 'mutbin' # 二进制变异
        else:
            self.mutFunc = 'mutpolyn' # 多项式变异
            if population.conordis == 0:
                self.recFunc = 'recsbx' # 模拟二进制交叉
            elif population.conordis == 1:
                self.recFunc = 'xovud' # 均匀交叉
        self.pc = 1 # 重组概率
        self.pm = 1 # 变异概率
        self.MAXSIZE = population.sizes # 非支配解集大小限制
        self.drawing = 1 # 绘图
        self.ax = None
        self.passTime = 0 # 记录用时
        
    def terminated(self, NDSet, population): # 判断是终止进化
        if self.currentGen < self.MAXGEN or NDSet.sizes > self.MAXSIZE:
            self.passTime += time.time() - self.timeSlot # 更新用时记录
            if self.drawing == 2:
                self.ax = ea.moeaplot(population.ObjV, False, self.ax, self.currentGen) # 绘制动态图
            self.timeSlot = time.time() # 更新时间戳
            self.currentGen += 1 # 进化代数+1
            return False
        else:
            return True
    
    def run(self):
        #==========================初始化配置===========================
        problem = self.problem
        population = self.population
        NIND = population.sizes
        MAXSIZE = self.MAXSIZE
        if MAXSIZE is None: # 检查MAXSIZE，默认取2倍的种群规模
            MAXSIZE = 2 * NIND
        aimFuc = problem.aimFuc # 获取目标函数地址
        #===========================准备进化============================
        self.timeSlot = time.time() # 开始计时
        if population.Chrom is None:
            population.initChrom(NIND) # 初始化种群染色体矩阵（内含解码，详见Population类的源码）
        population.ObjV, population.CV = aimFuc(population.Phen, population.CV) # 计算种群的目标函数值
        population.FitnV, NDSet = updateNDSet(population, problem.maxormins, MAXSIZE) # 计算适应度和得到全局非支配种群
        self.evalsNum = population.sizes # 记录评价次数
        #===========================开始进化============================
        self.currentGen = 0
        while self.terminated(NDSet, population) == False:
            uniChrom = np.unique(NDSet.Chrom, axis = 0)
            repRate = 1 - uniChrom.shape[0] / NDSet.sizes # 计算NDSet中的重复率
            # 选择基个体去进化形成子代
            offspring = population[ea.selecting(self.selFunc, population.FitnV, NIND)]
            # 对育种种群进行进化操作
            offspring.Chrom = ea.recombin(self.recFunc, offspring.Chrom, self.pc) #重组
            offspring.Chrom = ea.mutate(self.mutFunc, offspring.Encoding, offspring.Chrom, offspring.Field, self.pm) # 变异
            if population.Encoding != 'B' and population.Encoding != 'G' and repRate > 0.1:
                offspring.Chrom = ea.mutate('mutgau', offspring.Encoding, offspring.Chrom, offspring.Field, self.pm, False, 3) # 高斯变异，对标准差放大3倍。
            offspring.Phen = offspring.decoding() # 染色体解码
            offspring.ObjV, offspring.CV = aimFuc(offspring.Phen, offspring.CV) # 求进化后个体的目标函数值
            self.evalsNum += offspring.sizes # 更新评价次数
            # 父代种群和育种种群合并
            population = population + offspring
            population.FitnV, NDSet = updateNDSet(population, problem.maxormins, MAXSIZE, NDSet) # 计算合并种群的适应度及更新NDSet
            # 保留个体到下一代
            population = population[ea.selecting('dup', population.FitnV, NIND)] # 选择，保留NIND个个体
        NDSet = NDSet[np.where(np.all(NDSet.CV <= 0, 1))[0]] # 最后要彻底排除非可行解
        self.passTime += time.time() - self.timeSlot # 更新用时记录
        #=========================绘图及输出结果=========================
        if self.drawing != 0:
            ea.moeaplot(NDSet.ObjV, True)
        
        # 返回帕累托最优集以及执行时间
        return NDSet
    
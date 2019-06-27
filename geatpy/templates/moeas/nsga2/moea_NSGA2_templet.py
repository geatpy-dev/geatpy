# -*- coding: utf-8 -*-
import time
import numpy as np
import geatpy as ea # 导入geatpy库
from sys import path as paths
from os import path
paths.append(path.split(path.split(path.realpath(__file__))[0])[0])

class moea_NSGA2_templet(ea.Algorithm):
    
    """
moea_NSGA2_templet : class - 基于NSGA-II算法求解多目标优化问题的进化算法模板
    
算法描述:
    采用NSGA-II进行多目标优化，算法详见参考文献[1]。
    
模板使用注意:
    1.本模板调用的目标函数形如：[ObjV,CV] = aimFuc(Vars, CV), 
      其中Vars表示决策变量矩阵, CV为种群的违反约束程度矩阵(详见Geatpy数据结构),
    若不符合上述规范，则请修改算法模板或自定义新算法模板。

参考文献:
    [1] Deb K , Pratap A , Agarwal S , et al. A fast and elitist multiobjective 
    genetic algorithm: NSGA-II[J]. IEEE Transactions on Evolutionary 
    Computation, 2002, 6(2):0-197.
        
    """
    
    def __init__(self, problem, population):
        self.name = 'NSGA2'
        self.problem = problem
        self.population = population
        self.ndSort = ea.ndsortESS # 设置非支配排序算子
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
        self.drawing = 1 # 绘图
        self.ax = None
        self.passTime = 0 # 记录用时
    
    def calFitnV(self, population, NUM):
        """
        描述:
            计算种群个体的适应度
        算法:
            先对所需个体进行非支配排序分级，然后对所有个体进行拥挤度计算，最后先让处于较前帕累托分层的
            个体得到较高的适应度，对于同一层，则让拥挤距离大的个体得到较高的适应度
        输出参数:
            FitnV : array - 种群个体的适应度列向量
        """
        
        [levels, criLevel] = self.ndSort(self.problem.maxormins * population.ObjV, NUM, None, population.CV) # 对NUM个个体进行非支配分层
        dis = ea.crowdis(population.ObjV, levels) # 计算拥挤距离
        FitnV = np.array([np.argsort(np.lexsort(np.array([dis, -levels])), kind = 'mergesort')]).T # 计算适应度
        return FitnV
    
    def terminated(self, population): # 判断是否终止进化
        if self.currentGen < self.MAXGEN:
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
        self.ax = None # 存储上一桢动画
        population = self.population
        NIND = population.sizes
        #===========================准备进化============================
        self.timeSlot = time.time() # 开始计时
        if population.Chrom is None:
            population.initChrom(NIND) # 初始化种群染色体矩阵（内含解码，详见Population类的源码）
        population.ObjV, population.CV = self.problem.aimFuc(population.Phen, population.CV) # 计算种群的目标函数值
        self.evalsNum = population.sizes # 记录评价次数
        population.FitnV = self.calFitnV(population, NIND) # 计算适应度
        #===========================开始进化============================
        self.currentGen = 0
        while self.terminated(population) == False:
            # 选择基个体
            offspring = population[ea.selecting(self.selFunc, population.FitnV, NIND)]
            # 对基个体进行进化操作
            offspring.Chrom = ea.recombin(self.recFunc, offspring.Chrom, self.pc) #重组
            offspring.Chrom = ea.mutate(self.mutFunc, offspring.Encoding, offspring.Chrom, offspring.Field, self.pm) # 变异
            offspring.Phen = offspring.decoding() # 解码
            offspring.ObjV, offspring.CV = self.problem.aimFuc(offspring.Phen, offspring.CV) # 求进化后个体的目标函数值
            self.evalsNum += offspring.sizes # 更新评价次数
            # 合并
            population = population + offspring
            # 计算合并种群的适应度
            population.FitnV = self.calFitnV(population, NIND)
            # 选择个体保留到下一次进化
            population = population[ea.selecting('dup', population.FitnV, NIND)] # 调用低级选择算子dup进行基于适应度排序的选择，保留NIND个个体
        # 得到非支配种群
        [levels, criLevel] = self.ndSort(self.problem.maxormins * population.ObjV, NIND, 1, population.CV) # 非支配分层
        NDSet = population[np.where(levels == 1)[0]] # 只保留种群中的非支配个体，形成一个非支配种群
        NDSet = NDSet[np.where(np.all(NDSet.CV <= 0, 1))[0]] # 最后要彻底排除非可行解
        self.passTime += time.time() - self.timeSlot # 更新用时记录
        #=========================绘图及输出结果=========================
        if self.drawing != 0:
            ea.moeaplot(NDSet.ObjV, True)
        
        # 返回帕累托最优集
        return NDSet


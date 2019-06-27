# -*- coding: utf-8 -*-
import time
import numpy as np
import geatpy as ea # 导入geatpy库
from sys import path as paths
from os import path as path
paths.append(path.split(path.split(path.realpath(__file__))[0])[0])

class moea_NSGA3_templet(ea.Algorithm): # 继承Algorithm算法模板父类
    
    """
moea_NSGA3_templet : class - 基于NSGA-III算法求解多目标优化问题的进化算法模板
    
算法描述:
    采用NSGA-III进行多目标优化
    
模板使用注意:
    1.本模板调用的目标函数形如：[ObjV,CV] = aimFuc(Vars, CV), 
      其中Vars表示决策变量矩阵, CV为种群的可行度矩阵(详见Geatpy数据结构),
    若不符合上述规范，则请修改算法模板或自定义新算法模板。

参考文献:
    [1] Deb K , Jain H . An Evolutionary Many-Objective Optimization Algorithm 
    Using Reference-Point-Based Nondominated Sorting Approach, Part I: 
    Solving Problems With Box Constraints[J]. IEEE Transactions on 
    Evolutionary Computation, 2014, 18(4):577-601.
    
    """
    
    def __init__(self, problem, population):
        self.name = 'NSGA3'
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
        self.ax = None # 存储上一桢动画
        self.passTime = 0 # 记录用时
    
    def calFitnV(self, population, NUM, uniformPoint):
        """
        描述:
            计算种群个体的适应度
        算法:
            先对所需个体进行非支配排序分级，然后利用参考点关联选出所需要数量的点，最后结合帕累托分级以及是否被选择来计算适应度。    
        输出参数:
            FitnV : array - 种群个体的适应度列向量
        """
        
        [levels, criLevel] = self.ndSort(self.problem.maxormins * population.ObjV, NUM, None, population.CV) # 对NUM个个体进行非支配分层
        chooseFlag = ea.refselect(self.problem.maxormins * population.ObjV, levels, criLevel, NUM, uniformPoint, True) # 根据参考点选择个体(True表示使用伪随机数方法，可以提高速度，详见refselect帮助文档)
        FitnV = np.array([1 / (levels - chooseFlag + 1)]).T # 计算适应度
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
        population = self.population
        self.timeSlot = time.time() # 开始计时
        uniformPoint, NIND = ea.crtup(self.problem.M, population.sizes) # 生成在单位目标维度上均匀分布的参考点集
        if population.Chrom is None or population.sizes != NIND:
            population.initChrom(NIND) # 初始化种群染色体矩阵（内含解码，详见Population类的源码），此时种群规模将调整为uniformPoint点集的大小，initChrom函数会把种群规模给重置
        population.ObjV, population.CV = self.problem.aimFuc(population.Phen, population.CV) # 计算种群的目标函数值
        self.evalsNum = population.sizes # 记录评价次数
        #===========================开始进化============================
        self.currentGen = 0
        while self.terminated(population) == False:
            # 选择个体参与进化
            offspring = population[ea.selecting(self.selFunc, population.FitnV, NIND)]
            # 对基个体进行进化操作
            offspring.Chrom = ea.recombin(self.recFunc, offspring.Chrom, self.pc) #重组
            offspring.Chrom = ea.mutate(self.mutFunc, offspring.Encoding, offspring.Chrom, offspring.Field, self.pm) # 变异
            offspring.Phen = offspring.decoding() # 解码
            offspring.ObjV, offspring.CV = self.problem.aimFuc(offspring.Phen, offspring.CV) # 求进化后个体的目标函数值
            self.evalsNum += offspring.sizes # 更新评价次数
            # 合并
            population = population + offspring
            population.FitnV = self.calFitnV(population, NIND, uniformPoint) # 计算合并种群的适应度
            population = population[ea.selecting('dup', population.FitnV, NIND)] # 选择操作，保留NIND个个体
        # 得到非支配种群
        [levels, criLevel] = self.ndSort(self.problem.maxormins * population.ObjV, NIND, 1, population.CV) # 非支配分层
        NDSet = population[np.where(levels == 1)[0]] # 只保留种群中的非支配个体，形成一个非支配种群
        NDSet = NDSet[np.where(np.all(NDSet.CV <= 0, 1))[0]] # 最后要彻底排除非可行解
        self.passTime += time.time() - self.timeSlot # 更新用时记录
        # 绘图
        if self.drawing != 0:
            ea.moeaplot(NDSet.ObjV, True)
        
        # 返回帕累托最优集
        return NDSet

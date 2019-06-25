# -*- coding: utf-8 -*-
import time
import numpy as np
import geatpy as ea # 导入geatpy库
from sys import path as paths
from os import path as path
paths.append(path.split(path.split(path.realpath(__file__))[0])[0])

class moea_NSGA3_DE_templet(ea.Algorithm): # 继承Algorithm算法模板父类
    
    """
moea_NSGA3_DE_templet : class - 基于NSGA-III-DE算法求解多目标优化问题的进化算法模板
    
算法描述:
    采用NSGA-III-DE进行多目标优化，
    与NSGA-III不同的是，该算法把NSGA-III中的子代生成部分替换成DE/rand/1/bin。
    
模板使用注意:
    1.本模板调用的目标函数形如：[ObjV,CV] = aimFuc(Vars, CV), 
      其中Vars表示决策变量矩阵, CV为种群的可行度矩阵(详见Geatpy数据结构),
    若不符合上述规范，则请修改算法模板或自定义新算法模板。

参考文献:
    [1] Deb K , Jain H . An Evolutionary Many-Objective Optimization Algorithm 
    Using Reference-Point-Based Nondominated Sorting Approach, Part I: 
    Solving Problems With Box Constraints[J]. IEEE Transactions on 
    Evolutionary Computation, 2014, 18(4):577-601.
    
    [2] Tanabe R., Fukunaga A. (2014) Reevaluating Exponential Crossover in 
    Differential Evolution. In: Bartz-Beielstein T., Branke J., Filipič B., 
    Smith J. (eds) Parallel Problem Solving from Nature – PPSN XIII. PPSN 2014. 
    Lecture Notes in Computer Science, vol 8672. Springer, Cham
    
    """
    
    def __init__(self, problem, population):
        self.name = 'NSGA3-DE'
        self.problem = problem
        self.population = population
        self.ndSort = ea.ndsortESS # 设置非支配排序算子
        if population.Encoding == 'R' or population.Encoding == 'I':
            self.mutFunc = 'mutde' # 差分变异
            self.recFunc = 'xovbd' # 二项式分布交叉
        self.F = 0.5 # 差分变异缩放因子（可以设置为一个数也可以设置为一个列数与种群规模数目相等的列向量）
        self.pc = 0.2 # 交叉概率
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
        chooseFlag = ea.refselect(population.ObjV, levels, criLevel, NUM, uniformPoint, True) # 根据参考点选择个体(True表示使用伪随机数方法，可以提高速度，详见refselect帮助文档)
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
            population.shuffle() # 打乱个体顺序
            # 进行差分进化操作
            mutPop = population.copy()
            mutPop.Chrom = ea.mutate(self.mutFunc, mutPop.Encoding, mutPop.Chrom, mutPop.Field, self.F, 1) # 差分变异
            tempPop = population + mutPop # 当代种群个体与变异个体进行合并（为的是后面用于重组）
            experimentPop = population.copy() # 试验种群
            experimentPop.Chrom = ea.recombin(self.recFunc, tempPop.Chrom, self.pc, True) # 重组
            # 求进化后个体的目标函数值
            experimentPop.Phen = experimentPop.decoding() # 染色体解码
            experimentPop.ObjV, experimentPop.CV = self.problem.aimFuc(experimentPop.Phen, experimentPop.CV)
            self.evalsNum += experimentPop.sizes # 更新评价次数
            # 合并
            population = population + experimentPop
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

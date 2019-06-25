# -*- coding: utf-8 -*-
import time
import numpy as np
import geatpy as ea # 导入geatpy库
from sys import path as paths
from os import path
paths.append(path.split(path.split(path.realpath(__file__))[0])[0])

class soea_studGA_templet(ea.Algorithm):
    
    """
soea_studGA_templet.py - Stud GA templet(种马遗传算法模板)

算法描述:
    本模板实现的是种马遗传算法。算法流程详见参考文献[1]。
    为了使算法也能很好地求解约束优化问题，本算法模板稍作修改，增添“遗忘策略”，
    当某一代没有可行个体时，让进化记录器忽略这一代，不对这一代的个体进行记录。
    
模板使用注意：
    1.本模板调用的目标函数形如：[ObjV,CV] = aimFuc(Vars, CV), 
      其中Vars表示决策变量矩阵, CV为种群的违反约束程度矩阵(详见Geatpy数据结构),
    若不符合上述规范，则请修改算法模板或自定义新算法模板。

参考文献:
    [1] Khatib W , Fleming P J . The stud GA: A mini revolution?[C]// International 
    Conference on Parallel Problem Solving from Nature. Springer, Berlin, Heidelberg, 1998.
    
"""
    
    def __init__(self, problem, population):
        self.name = 'studGA'
        self.problem = problem
        self.population = population
        self.selFunc = 'tour' # 锦标赛选择算子
        if population.Encoding == 'B' or population.Encoding == 'G':
            self.recFunc = 'xovdp' # 两点交叉
            self.mutFunc = 'mutbin' # 二进制变异
        elif population.Encoding == 'R' or population.Encoding == 'I':
            self.recFunc = 'xovdp' # 两点交叉
            if population.conordis == 0:
                self.mutFunc = 'mutbga' # 实数值变异
            elif population.conordis == 1:
                self.mutFunc = 'mutuni' # 均匀变异
        elif population.Encoding == 'P':
            self.recFunc = 'xovpmx' # 部分匹配交叉
            self.mutFunc = 'mutpp' # 排列编码染色体变异
        self.pc = 1 # 重组概率
        self.pm = 1 # 变异概率
        self.drawing = 1 # 绘图
        self.ax = None # 存储上一桢动画
        self.passTime = 0 # 记录用时
    
    def stat(self, population): # 分析记录
        # 进行进化记录
        feasible = np.where(np.all(population.CV <= 0, 1))[0] # 找到可行解的下标
        if len(feasible) > 0:
            tempPop = population[feasible]
            bestIdx = np.argmax(tempPop.FitnV) # 获取最优个体的下标 
            self.obj_trace[self.currentGen,0] = np.sum(tempPop.ObjV) / tempPop.sizes # 记录种群个体平均目标函数值
            self.obj_trace[self.currentGen,1] = tempPop.ObjV[bestIdx] # 记录当代目标函数的最优值
            self.var_trace[self.currentGen,:] = tempPop.Phen[bestIdx, :] # 记录当代最优的决策变量值
            self.forgetCount = 0 # “遗忘策略”计数器清零
            self.passTime += time.time() - self.timeSlot # 更新用时记录
            if self.drawing == 2:
                self.ax = ea.soeaplot(self.obj_trace[:,[1]],'种群最优个体目标函数值', False, self.ax, self.currentGen) # 绘制动态图
            self.timeSlot = time.time() # 更新时间戳
        else:
            self.currentGen -= 1 # 忽略这一代
            self.forgetCount += 1 # “遗忘策略”计数器加1
    
    def terminated(self, population): # 判断是否终止进化
        if self.currentGen < self.MAXGEN or self.forgetCount >= 10 * self.MAXGEN:
            self.stat(population) # 分析记录当代种群的数据
            self.currentGen += 1 # 进化代数+1
            return False
        else:
            return True
    
    def run(self):
        #==========================初始化配置===========================
        population = self.population
        NIND = population.sizes
        NVAR = self.problem.Dim # 得到决策变量的个数
        self.obj_trace = (np.zeros((self.MAXGEN, 2)) * np.nan) # 定义目标函数值记录器，初始值为nan
        self.var_trace = (np.zeros((self.MAXGEN, NVAR)) * np.nan) # 定义变量记录器，记录决策变量值，初始值为nan
        self.forgetCount = 0 # “遗忘策略”计数器，用于记录连续出现最优个体不是可行个体的代数
        #===========================准备进化============================
        self.timeSlot = time.time() # 开始计时
        if population.Chrom is None:
            population.initChrom(NIND) # 初始化种群染色体矩阵（内含染色体解码，详见Population类的源码）
        population.ObjV, population.CV = self.problem.aimFuc(population.Phen, population.CV) # 计算种群的目标函数值
        population.FitnV = ea.scaling(self.problem.maxormins * population.ObjV, population.CV) # 计算适应度
        self.evalsNum = population.sizes # 记录评价次数
        #===========================开始进化============================
        self.currentGen = 0
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
            population.ObjV, population.CV = self.problem.aimFuc(population.Phen, population.CV)
            self.evalsNum += population.sizes # 更新评价次数
            population.FitnV = ea.scaling(self.problem.maxormins * population.ObjV, population.CV) # 计算适应度
        # 处理进化记录器
        delIdx = np.where(np.isnan(self.obj_trace))[0]
        self.obj_trace = np.delete(self.obj_trace, delIdx, 0)
        self.var_trace = np.delete(self.var_trace, delIdx, 0)
        if self.obj_trace.shape[0] == 0:
            raise RuntimeError('error: No feasible solution. (有效进化代数为0，没找到可行解。)')
        self.passTime += time.time() - self.timeSlot # 更新用时记录
        # 绘图
        if self.drawing != 0:
            ea.trcplot(self.obj_trace, [['种群个体平均目标函数值', '种群最优个体目标函数值']])
        # 返回最后一代种群、进化记录器、变量记录器以及执行时间
        return [population, self.obj_trace, self.var_trace]
    
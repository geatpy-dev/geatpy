# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea
import warnings
import time

class Algorithm:
    
    """
Algorithm : class - 算法模板顶级父类

描述:
    算法设置类是用来存储与算法运行参数设置相关信息的一个类。

属性:
    name            : str      - 算法名称（可以自由设置名称）。
    
    problem         : class <Problem> - 问题类的对象。
    
    MAXGEN          : int      - 最大进化代数。
    
    currentGen      : int      - 当前进化的代数。
    
    MAXTIME         : float    - 时间限制（单位：秒）。
    
    timeSlot        : float    - 时间戳（单位：秒）。
    
    passTime        : float    - 已用时间（单位：秒）。
    
    MAXEVALS        : int      - 最大评价次数。
    
    evalsNum        : int      - 当前评价次数。
    
    MAXSIZE         : int      - 最优解的最大数目。
        
    population      : class <Population> - 种群对象。
    
    drawing         : int      - 绘图方式的参数，
                                 0表示不绘图，
                                 1表示绘制结果图，
                                 2表示实时绘制目标空间动态图，
                                 3表示实时绘制决策空间动态图。

函数:
    terminated()    : 计算是否需要终止进化，具体功能需要在继承类即算法模板中实现。
    
    run()           : 执行函数，需要在继承类即算法模板中实现。
    
    check()         : 用于检查种群对象的ObjV和CV的数据是否有误。
    
    call_aimFunc()  : 用于调用问题类中的aimFunc()进行计算ObjV和CV(若有约束)。
    
"""

    def __init__(self):
        self.name = 'Algorithm'
        self.problem = None
        self.MAXGEN = None
        self.currentGen = None
        self.MAXTIME = None
        self.timeSlot = None
        self.passTime = None
        self.MAXEVALS = None
        self.evalsNum = None
        self.MAXSIZE = None
        self.population = None
        self.drawing = None
    
    def terminated(self):
        pass
    
    def run(self):
        pass
    
    def check(self, pop):
        
        """
        用于检查种群对象的ObjV和CV的数据是否有误。
        
        """
        
        # 检测数据非法值
        if np.any(np.isnan(pop.ObjV)):
            warnings.warn("Warning: Some elements of ObjV are NAN, please check the calculation of ObjV.(ObjV的部分元素为NAN，请检查目标函数的计算。)", RuntimeWarning)
        elif np.any(np.isinf(pop.ObjV)):
            warnings.warn("Warning: Some elements of ObjV are Inf, please check the calculation of ObjV.(ObjV的部分元素为Inf，请检查目标函数的计算。)", RuntimeWarning)
        if pop.CV is not None:
            if np.any(np.isnan(pop.CV)):
                warnings.warn("Warning: Some elements of CV are NAN, please check the calculation of CV.(CV的部分元素为NAN，请检查CV的计算。)", RuntimeWarning)
            elif np.any(np.isinf(pop.CV)):
                warnings.warn("Warning: Some elements of CV are Inf, please check the calculation of CV.(CV的部分元素为Inf，请检查CV的计算。)", RuntimeWarning)
    
    def call_aimFunc(self, pop):
        
        """
        使用注意:
        本函数调用的目标函数形如：aimFunc(pop), (在自定义问题类中实现)。
        其中pop为种群类的对象，代表一个种群，
        pop对象的Phen属性（即种群染色体的表现型）等价于种群所有个体的决策变量组成的矩阵，
        该函数根据该Phen计算得到种群所有个体的目标函数值组成的矩阵，并将其赋值给pop对象的ObjV属性。
        若有约束条件，则在计算违反约束程度矩阵CV后赋值给pop对象的CV属性（详见Geatpy数据结构）。
        该函数不返回任何的返回值，求得的目标函数值保存在种群对象的ObjV属性中，
                              违反约束程度矩阵保存在种群对象的CV属性中。
        例如：population为一个种群对象，则调用call_aimFunc(population)即可完成目标函数值的计算。
             之后可通过population.ObjV得到求得的目标函数值，population.CV得到违反约束程度矩阵。
        若不符合上述规范，则请修改算法模板或自定义新算法模板。
        
        """
        
        pop.Phen = pop.decoding() # 染色体解码
        if self.problem is None:
            raise RuntimeError('error: problem has not been initialized. (算法模板中的问题对象未被初始化。)')
        self.problem.aimFunc(pop) # 调用问题类的aimFunc()
        self.evalsNum = self.evalsNum + pop.sizes if self.evalsNum is not None else pop.sizes # 更新评价次数
        if type(pop.ObjV) != np.ndarray or pop.ObjV.ndim != 2 or pop.ObjV.shape[0] != pop.sizes or pop.ObjV.shape[1] != self.problem.M:
            raise RuntimeError('error: ObjV is illegal. (目标函数值矩阵ObjV的数据格式不合法，请检查目标函数的计算。)')
        if pop.CV is not None:
            if type(pop.CV) != np.ndarray or pop.CV.ndim != 2 or pop.CV.shape[0] != pop.sizes:
                raise RuntimeError('error: CV is illegal. (违反约束程度矩阵CV的数据格式不合法，请检查CV的计算。)')
        
class MoeaAlgorithm(Algorithm): # 多目标优化算法模板父类
    
    """
    描述:
        此为多目标进化优化算法模板的父类，所有多目标优化算法模板均继承自该父类。
        为了使算法也能很好地求解约束优化问题，本算法模板稍作修改，增添“遗忘策略”，
        当某一代没有可行个体时，让进化记录器忽略这一代，不对这一代的个体进行记录，但不影响进化。
        
    """
    
    def __init__(self, problem, population): # 构造方法，这里只初始化静态参数以及对动态参数进行定义
        Algorithm.__init__(self) # 先调用父类构造方法
        self.problem = problem
        self.population = population
        self.drawing = 1 # 绘图
        self.ax = None # 用于存储动态图
        self.forgetCount = None # “遗忘策略”计数器，用于记录连续若干代出现种群所有个体都不是可行个体的次数
        self.maxForgetCount = None # “遗忘策略”计数器最大上限值
        self.pop_trace = None # 种群记录器
    
    def initialization(self):
        
        """
        描述: 该函数用于在进化前对算法模板的参数进行初始化操作。
        该函数需要在执行算法模板的run()方法的一开始被调用，同时开始计时，
        以确保所有这些参数能够被正确初始化。
        
        """
        self.ax = None # 重置ax
        self.passTime = 0 # 初始化计时器
        self.forgetCount = 0 # 初始化“遗忘策略”计数器
        self.maxForgetCount = 1000 # 初始化“遗忘策略”计数器最大上限值，当超过这个上限时将终止进化
        self.pop_trace = [] # 初始化种群记录器
        self.currentGen = 0 # 设置初始为第0代
        self.evalsNum = 0 # 设置评价次数为0
        self.timeSlot = time.time() # 开始计时
    
    def stat(self, population): # 分析记录，更新进化记录器，population为传入的种群对象
        feasible = np.where(np.all(population.CV <= 0, 1))[0] if population.CV is not None else np.arange(population.sizes) # 找到可行解个体的下标
        if len(feasible) > 0:
            tempPop = population[feasible] # 获取可行解个体
            self.pop_trace.append(tempPop) # 添加记录（只添加可行解个体到种群记录器中）
            self.forgetCount = 0 # “遗忘策略”计数器清零
            self.passTime += time.time() - self.timeSlot # 更新用时记录
            if self.drawing == 2:
                # 绘制目标空间动态图
                self.ax = ea.moeaplot(tempPop.ObjV, 'objective values', False, self.ax, self.currentGen, gridFlag = True)
            elif self.drawing == 3:
                # 绘制决策空间动态图
                self.ax = ea.varplot(tempPop.Phen, 'decision variables', False, self.ax, self.currentGen, gridFlag = False)
            self.timeSlot = time.time() # 更新时间戳
        else:
            self.currentGen -= 1 # 忽略这一代
            self.forgetCount += 1 # “遗忘策略”计数器加1
        
    def terminated(self, population):
        
        """
        描述:
            该函数用于判断是否应该终止进化，population为传入的种群。
            
        """
        
        self.check(population) # 检查种群对象的关键属性是否有误
        self.stat(population) # 进行统计分析，更新进化记录器
        # 判断是否终止进化，由于代数是从0数起，因此在比较currentGen和MAXGEN时需要对currentGen加1
        if self.currentGen + 1 >= self.MAXGEN or self.forgetCount >= self.maxForgetCount:
            return True
        else:
            self.currentGen += 1 # 进化代数+1
            return False
    
    def finishing(self, population):
        
        """
        进化完成后调用的函数。
        
        """
        
        # 得到非支配种群
        [levels, criLevel] = ea.ndsortDED(population.ObjV, None, 1, population.CV, self.problem.maxormins) # 非支配分层
        NDSet = population[np.where(levels == 1)[0]] # 只保留种群中的非支配个体，形成一个非支配种群
        if NDSet.CV is not None: # CV不为None说明有设置约束条件
            NDSet = NDSet[np.where(np.all(NDSet.CV <= 0, 1))[0]] # 最后要彻底排除非可行解
        self.passTime += time.time() - self.timeSlot # 更新用时记录
        if NDSet.sizes == 0:
            raise RuntimeError('error: No feasible solution. (没找到可行解。)')
        # 绘图
        if self.drawing != 0:
            if NDSet.ObjV.shape[1] == 2 or NDSet.ObjV.shape[1] == 3:
                ea.moeaplot(NDSet.ObjV, 'Pareto Front', saveFlag = True, gridFlag = True)
            else:
                ea.moeaplot(NDSet.ObjV, 'Value Path', saveFlag = True, gridFlag = False)
        # 返回帕累托最优集
        return NDSet

class SoeaAlgorithm(Algorithm): # 单目标优化算法模板父类
    
    """
    描述:
        此为单目标进化优化算法模板的父类，所有单目标优化算法模板均继承自该父类。
        为了使算法也能很好地求解约束优化问题，本算法模板稍作修改，增添“遗忘策略”，
        当某一代没有可行个体时，让进化记录器忽略这一代，不对这一代的个体进行记录，但不影响进化。
    """
    
    def __init__(self, problem, population): # 构造方法，这里只初始化静态参数以及对动态参数进行定义
        Algorithm.__init__(self) # 先调用父类构造方法
        self.problem = problem
        self.population = population
        self.drawing = 1 # 绘图
        self.forgetCount = None # “遗忘策略”计数器，用于记录连续若干代出现种群所有个体都不是可行个体的次数
        self.maxForgetCount = 1000 # “遗忘策略”计数器最大上限值，当超过这个上限时将终止进化
        self.trappedCount = 0 # “进化停滞”计数器
        self.trappedValue = 0 # 进化算法陷入停滞的判断阈值，当abs(最优目标函数值-上代的目标函数值) < trappedValue时，对trappedCount加1
        self.maxTrappedCount = 1000 # 进化停滞计数器最大上限值，当超过这个上限时将终止进化
        self.preObjV = np.nan # “前代最优目标函数值记录器”，用于记录上一代的最优目标函数值
        self.ax = None # 存储上一桢动画
    
    def initialization(self):
        
        """
        描述: 该函数用于在进化前对算法模板的一些动态参数进行初始化操作
        该函数需要在执行算法模板的run()方法的一开始被调用，同时开始计时，
        以确保所有这些参数能够被正确初始化。
        
        """
        
        self.ax = None # 重置ax
        self.passTime = 0 # 记录用时
        self.forgetCount = 0 # “遗忘策略”计数器，用于记录连续若干代出现种群所有个体都不是可行个体的次数
        self.preObjV = np.nan # 重置“前代最优目标函数值记录器”
        self.trappedCount = 0 # 重置“进化停滞”计数器
        self.obj_trace = np.zeros((self.MAXGEN, 2)) * np.nan # 定义目标函数值记录器，初始值为nan
        self.var_trace = np.zeros((self.MAXGEN, self.problem.Dim)) * np.nan # 定义变量记录器，记录决策变量值，初始值为nan
        self.currentGen = 0 # 设置初始为第0代
        self.evalsNum = 0 # 设置评价次数为0
        self.timeSlot = time.time() # 开始计时

    def stat(self, population): # 分析记录，更新进化记录器，population为传入的种群对象
        # 进行进化记录
        feasible = np.where(np.all(population.CV <= 0, 1))[0] if population.CV is not None else np.arange(population.sizes) # 找到可行解个体的下标
        if len(feasible) > 0:
            tempPop = population[feasible]
            bestIdx = np.argmax(tempPop.FitnV) # 获取最优个体的下标
            self.obj_trace[self.currentGen, 0] = np.sum(tempPop.ObjV) / tempPop.sizes # 记录种群个体平均目标函数值
            self.obj_trace[self.currentGen, 1] = tempPop.ObjV[bestIdx] # 记录当代目标函数的最优值
            self.var_trace[self.currentGen, :] = tempPop.Phen[bestIdx, :] # 记录当代最优的决策变量值
            self.forgetCount = 0 # “遗忘策略”计数器清零
            if np.abs(self.preObjV - self.obj_trace[self.currentGen, 1]) < self.trappedValue:
                self.trappedCount += 1
            else:
                self.trappedCount = 0 # 重置进化停滞计数器
            self.passTime += time.time() - self.timeSlot # 更新用时记录
            if self.drawing == 2:
                self.ax = ea.soeaplot(self.obj_trace[:,[1]], Label = 'Objective Value' , saveFlag = False, ax = self.ax, gen = self.currentGen, gridFlag = False) # 绘制动态图
            elif self.drawing == 3:
                self.ax = ea.varplot(tempPop.Phen, Label = 'decision variables', saveFlag = False, ax = self.ax, gen = self.currentGen, gridFlag = False)
            self.timeSlot = time.time() # 更新时间戳
        else:
            self.currentGen -= 1 # 忽略这一代
            self.forgetCount += 1 # “遗忘策略”计数器加1
    
    def terminated(self, population):
        
        """
        描述:
            该函数用于判断是否应该终止进化，population为传入的种群。
            
        """
        
        self.check(population) # 检查种群对象的关键属性是否有误
        self.stat(population) # 分析记录当代种群的数据
        # 判断是否终止进化，由于代数是从0数起，因此在比较currentGen和MAXGEN时需要对currentGen加1
        if self.currentGen + 1 >= self.MAXGEN or self.forgetCount >= self.maxForgetCount or self.trappedCount >= self.maxTrappedCount:
            return True
        else:
            self.preObjV = self.obj_trace[self.currentGen, 1] # 更新“前代最优目标函数值记录器”
            self.currentGen += 1 # 进化代数+1
            return False

    def finishing(self, population):
        
        """
        进化完成后调用的函数。
        
        """
        
        # 处理进化记录器
        delIdx = np.where(np.isnan(self.obj_trace))[0]
        self.obj_trace = np.delete(self.obj_trace, delIdx, 0)
        self.var_trace = np.delete(self.var_trace, delIdx, 0)
        if self.obj_trace.shape[0] == 0:
            raise RuntimeError('error: No feasible solution. (有效进化代数为0，没找到可行解。)')
        self.passTime += time.time() - self.timeSlot # 更新用时记录
        # 绘图
        if self.drawing != 0:
            ea.trcplot(self.obj_trace, [['种群个体平均目标函数值', '种群最优个体目标函数值']], xlabels = [['Number of Generation']], ylabels = [['Value']], gridFlags = [[False]])
        # 返回最后一代种群、进化记录器、变量记录器以及执行时间
        return [population, self.obj_trace, self.var_trace]
    
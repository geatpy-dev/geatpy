# -*- coding: utf-8 -*-
import time
import warnings
import numpy as np
import geatpy as ea


class Algorithm:
    """
Algorithm : class - 算法类的顶级父类

描述:
    算法设置类是用来存储与算法运行参数设置相关信息的一个类。

属性:
    name            : str      - 算法名称（可以自由设置名称）。
    
    problem         : class <Problem> - 问题类的对象。

    population      : class <Population> - 种群对象。
    
    MAXGEN          : int      - 最大进化代数。
    
    currentGen      : int      - 当前进化的代数。
    
    MAXTIME         : float    - 时间限制（单位：秒）。
    
    timeSlot        : float    - 时间戳（单位：秒）。
    
    passTime        : float    - 已用时间（单位：秒）。
    
    MAXEVALS        : int      - 最大评价次数。
    
    evalsNum        : int      - 当前评价次数。
    
    MAXSIZE         : int      - 最优个体的最大数目。

    outFunc         : function - 提供给用户自定义的函数，它在每发生一次进化时被调用。如果没有定义，则默认为None。
    
    logTras         : int      - Tras即周期的意思，该参数用于设置在进化过程中每多少代记录一次日志信息。
                                 设置为0表示不记录日志信息。
                                 注：此时假如设置了“每10代记录一次日志”而导致最后一代没有被记录，
                                     则会补充记录最后一代的信息，除非找不到可行解。

    log             : Dict     - 日志记录。其中包含2个基本的键：'gen'和'eval'，其他键的定义由该算法类的子类实现。
                                 'gen'的键值为一个list列表，用于存储日志记录中的每一条记录对应第几代种群。
                                 'eval'的键值为一个list列表，用于存储进化算法的评价次数。
                                 注：若设置了logTras为0，则不会记录日志，此时log会被设置为None。
    
    verbose         : bool     - 表示是否在输入输出流中打印输出日志信息。

    stopMsg         : str      - 记录进化终止原因的字符串。

    dirName         : str      - 用于指明文件保存的路径。用于把绘图文件保存在此目录下。
                                 当缺省或为None时，默认dirName=''，此时如果绘制图片，图片会被保存在执行文件的所在目录下。

函数:
    __init__()       : 构造函数，定义一些属性，并初始化一些静态参数。

    initialization() : 在进化前对算法类的一些动态参数进行初始化操作，具体功能需要在继承类中实现。
    
    run()            : 执行函数，具体功能需要在继承类中实现。

    logging()        : 用于在进化过程中记录日志，具体功能需要在继承类中实现。

    stat()           : 用于分析当代种群的信息，具体功能需要在继承类中实现。

    terminated()     : 计算是否需要终止进化，具体功能需要在继承类中实现。

    finishing ()     : 进化完成后调用的函数，具体功能需要在继承类中实现。

    check()          : 用于检查种群对象的ObjV和CV的数据是否有误。

    call_aimFunc()   : 用于调用问题类中的aimFunc()或evalVars()进行计算ObjV和CV(若有约束)。

    display()        : 用于在进化过程中进行一些输出，需要依赖属性verbose和log属性。

"""

    def __init__(self,
                 problem,
                 population,
                 MAXGEN = None,
                 MAXTIME = None,
                 MAXEVALS = None,
                 MAXSIZE = None,
                 logTras = None,
                 verbose = None,
                 outFunc = None,
                 dirName = None,
                 **kwargs):

        """
        描述: 
            构造函数。
        
        """

        if dirName is None:
            dirName = ''
        if dirName != '':
            dirName += '/'

        # 静态属性
        self.name = 'Algorithm'
        self.problem = problem
        self.population = population
        self.MAXGEN = MAXGEN
        self.MAXTIME = MAXTIME
        self.MAXEVALS = MAXEVALS
        self.MAXSIZE = MAXSIZE
        self.logTras = 1 if logTras is None else logTras
        self.verbose = True if verbose is None else verbose
        self.outFunc = outFunc
        self.dirName = dirName
        # 动态属性
        self.currentGen = None
        self.timeSlot = None
        self.passTime = None
        self.evalsNum = None
        self.log = None
        self.stopMsg = ''
        # 初始化种群的译码矩阵
        if isinstance(population, ea.Population):
            if population.ChromNum == 1:  # 如果是单染色体种群
                if population.Field is None:
                    self.population = ea.Population(Encoding=population.Encoding,
                                                    Field=(problem.varTypes, problem.ranges, problem.borders),
                                                    NIND=population.sizes)
            else:  # 如果是多染色体种群
                if population.Fields is None:
                    self.population = ea.PsyPopulation(Encodings=population.Encodings,
                                                       Fields=(problem.varTypes, problem.ranges, problem.borders),
                                                       NIND=population.sizes,
                                                       EncoIdxs=population.EncoIdxs)

    def initialization(self):
        pass

    def run(self, pop):
        pass

    def logging(self, pop):
        pass

    def stat(self, pop):
        pass

    def terminated(self, pop):
        pass

    def finishing(self, pop):
        pass

    def check(self, pop):

        """
        描述:
            用于检查种群对象的ObjV和CV的数据是否有误。

        输入参数:
            pop : class <Population> - 种群对象。

        输出参数:
            无输出参数。

        """

        # 检测数据非法值
        if np.any(np.isnan(pop.ObjV)):
            warnings.warn(
                "Warning: Some elements of ObjV are NAN, please check the calculation of ObjV.("
                "ObjV的部分元素为NAN，请检查目标函数的计算。)",
                RuntimeWarning)
        elif np.any(np.isinf(pop.ObjV)):
            warnings.warn(
                "Warning: Some elements of ObjV are Inf, please check the calculation of ObjV.("
                "ObjV的部分元素为Inf，请检查目标函数的计算。)",
                RuntimeWarning)
        if pop.CV is not None:
            if np.any(np.isnan(pop.CV)):
                warnings.warn(
                    "Warning: Some elements of CV are NAN, please check the calculation of CV.(CV的部分元素为NAN，请检查CV的计算。)",
                    RuntimeWarning)
            elif np.any(np.isinf(pop.CV)):
                warnings.warn(
                    "Warning: Some elements of CV are Inf, please check the calculation of CV.(CV的部分元素为Inf，请检查CV的计算。)",
                    RuntimeWarning)

    def call_aimFunc(self, pop):

        """
        描述: 调用问题类的aimFunc()或evalVars()完成种群目标函数值和违反约束程度的计算。

        例如：population为一个种群对象，则调用call_aimFunc(population)即可完成目标函数值的计算。
             之后可通过population.ObjV得到求得的目标函数值，population.CV得到违反约束程度矩阵。

        输入参数:
            pop : class <Population> - 种群对象。

        输出参数:
            无输出参数。

        """

        pop.Phen = pop.decoding()  # 染色体解码
        if self.problem is None:
            raise RuntimeError('error: problem has not been initialized. (算法类中的问题对象未被初始化。)')
        self.problem.evaluation(pop)  # 调用问题类的evaluation()
        self.evalsNum = self.evalsNum + pop.sizes if self.evalsNum is not None else pop.sizes  # 更新评价次数
        # 格式检查
        if not isinstance(pop.ObjV, np.ndarray) or pop.ObjV.ndim != 2 or pop.ObjV.shape[0] != pop.sizes or \
                pop.ObjV.shape[1] != self.problem.M:
            raise RuntimeError('error: ObjV is illegal. (目标函数值矩阵ObjV的数据格式不合法，请检查目标函数的计算。)')
        if pop.CV is not None:
            if not isinstance(pop.CV, np.ndarray) or pop.CV.ndim != 2 or pop.CV.shape[0] != pop.sizes:
                raise RuntimeError('error: CV is illegal. (违反约束程度矩阵CV的数据格式不合法，请检查CV的计算。)')

    def display(self):

        """
        描述:
            该函数打印日志log中每个键值的最后一条数据。假如log中只有一条数据或没有数据，则会打印表头。
            该函数将会在子类中被覆盖，以便进行更多其他的输出展示。

        """

        self.passTime += time.time() - self.timeSlot  # 更新用时记录，不计算display()的耗时
        headers = []
        widths = []
        values = []
        for key in self.log.keys():
            # 设置单元格宽度
            if key == 'gen':
                if self.MAXGEN is None:
                    width = 5
                else:
                    width = max(3, len(str(self.MAXGEN - 1)))  # 因为字符串'gen'长度为3，所以最小要设置长度为3
            elif key == 'eval':
                width = 8  # 因为字符串'eval'长度为4，所以最小要设置长度为4
            else:
                width = 13  # 预留13位显示长度，若数值过大，表格将无法对齐，此时若要让表格对齐，需要自定义算法类重写该函数
            headers.append(key)
            widths.append(width)
            value = self.log[key][-1] if len(self.log[key]) != 0 else "-"
            if isinstance(value, float):
                values.append("%.5E" % value)  # 格式化浮点数，输出时只保留至小数点后5位
            else:
                values.append(value)
        if len(self.log['gen']) == 1:  # 打印表头
            header_regex = '|'.join(['{}'] * len(headers))
            header_str = header_regex.format(*[str(key).center(width) for key, width in zip(headers, widths)])
            print("=" * len(header_str))
            print(header_str)
            print("-" * len(header_str))
        if len(self.log['gen']) != 0:  # 打印表格最后一行
            value_regex = '|'.join(['{}'] * len(values))
            value_str = value_regex.format(*[str(value).center(width) for value, width in zip(values, widths)])
            print(value_str)
        self.timeSlot = time.time()  # 更新时间戳


class MoeaAlgorithm(Algorithm):  # 多目标优化算法类的父类

    """
    描述:
        此为多目标进化优化算法类的父类，所有多目标优化算法类均继承自该父类。

    对比于父类该类新增的变量和函数:

        drawing        : int - 绘图方式的参数，
                               0表示不绘图；
                               1表示绘制最终结果图；
                               2表示实时绘制目标空间动态图；
                               3表示实时绘制决策空间动态图。

        draw()         : 绘图函数。

    """

    def __init__(self,
                 problem,
                 population,
                 MAXGEN = None,
                 MAXTIME = None,
                 MAXEVALS = None,
                 MAXSIZE = None,
                 logTras = None,
                 verbose = None,
                 outFunc = None,
                 drawing = None,
                 dirName = None,
                 **kwargs):

        """
        描述: 
            在该构造函数里只初始化静态参数以及对动态参数进行定义。
        
        """

        super().__init__(problem, population, MAXGEN, MAXTIME, MAXEVALS, MAXSIZE, logTras, verbose, outFunc, dirName)  # 先调用父类构造函数
        self.drawing = 1 if drawing is None else drawing
        # 以下为用户不需要设置的属性
        self.plotter = None  # 存储绘图对象

    def __str__(self):
        info = {}
        info['Algorithm Name'] = self.name
        info['Algorithm MAXGEN'] = self.MAXGEN
        info['Algorithm MAXTIME'] = self.MAXTIME
        info['Algorithm MAXEVALS'] = self.MAXEVALS
        return str(info)

    def initialization(self):

        """
        描述: 
            该函数用于在进化前对算法类的一些动态参数进行初始化操作。
            该函数需要在执行算法类的run()函数的一开始被调用，同时开始计时，
            以确保所有这些参数能够被正确初始化。
        
        """

        self.ax = None  # 初始化ax
        self.passTime = 0  # 初始化passTime
        self.log = None  # 初始化log
        self.currentGen = 0  # 初始为第0代
        self.evalsNum = 0  # 初始化评价次数为0
        self.log = {'gen': [], 'eval': []} if self.logTras != 0 else None  # 初始化log
        self.timeSlot = time.time()  # 开始计时

    def logging(self, pop):

        """
        描述:
            用于在进化过程中记录日志。该函数在stat()函数里面被调用。
            如果需要在日志中记录其他数据，需要在自定义算法类中重写该函数。

        输入参数:
            pop : class <Population> - 种群对象。

        输出参数:
            无输出参数。

        """

        self.passTime += time.time() - self.timeSlot  # 更新用时记录，不计算logging的耗时
        if len(self.log['gen']) == 0:  # 初始化log的各个键值
            self.log['gd'] = []
            self.log['igd'] = []
            self.log['hv'] = []
            self.log['spacing'] = []
        self.log['gen'].append(self.currentGen)
        self.log['eval'].append(self.evalsNum)  # 记录评价次数
        [levels, _] = ea.ndsortDED(pop.ObjV, needLevel=1, CV=pop.CV, maxormins=self.problem.maxormins)  # 非支配分层
        NDSet = pop[np.where(levels == 1)[0]]  # 只保留种群中的非支配个体，形成一个非支配种群
        if self.problem.ReferObjV is not None:
            self.log['gd'].append(ea.indicator.GD(NDSet.ObjV, self.problem.ReferObjV))  # 计算GD指标
            self.log['igd'].append(ea.indicator.IGD(NDSet.ObjV, self.problem.ReferObjV))  # 计算IGD指标
            self.log['hv'].append(ea.indicator.HV(NDSet.ObjV, self.problem.ReferObjV))  # 计算HV指标
        else:
            self.log['gd'].append(None)
            self.log['igd'].append(None)
            self.log['hv'].append(ea.indicator.HV(NDSet.ObjV))  # 计算HV指标
        self.log['spacing'].append(ea.indicator.Spacing(NDSet.ObjV))  # 计算Spacing指标
        self.timeSlot = time.time()  # 更新时间戳

    def draw(self, pop, EndFlag=False):

        """
        描述:
            该函数用于在进化过程中进行绘图。该函数在stat()以及finishing函数里面被调用。

        输入参数:
            pop     : class <Population> - 种群对象。
            
            EndFlag : bool - 表示是否是最后一次调用该函数。

        输出参数:
            无输出参数。

        """

        if not EndFlag:
            self.passTime += time.time() - self.timeSlot  # 更新用时记录，不计算画图的耗时
            # 绘制动画
            if self.drawing == 2:
                # 绘制目标空间动态图
                if pop.ObjV.shape[1] == 2 or pop.ObjV.shape[1] == 3:
                    if self.plotter is None:
                        self.plotter = ea.PointScatter(self.problem.M, grid=True, legend=True, title='Pareto Front Plot')
                    self.plotter.refresh()
                    self.plotter.add(pop.ObjV, color='red', label='MOEA PF at ' + str(self.currentGen) + ' Generation')
                else:
                    if self.plotter is None:
                        self.plotter = ea.ParCoordPlotter(self.problem.M, grid=True, legend=True, title='Parallel Coordinate Plot')
                    self.plotter.refresh()
                    self.plotter.add(pop.ObjV, color='red', label='MOEA Objective Value at ' + str(self.currentGen) + ' Generation')
                self.plotter.draw()
            elif self.drawing == 3:
                # 绘制决策空间动态图
                if self.plotter is None:
                    self.plotter = ea.ParCoordPlotter(self.problem.Dim, grid=True, legend=True, title='Variables Value Plot')
                self.plotter.refresh()
                self.plotter.add(pop.Phen, marker='o', color='blue', label='Variables Value at ' + str(self.currentGen) + ' Generation')
                self.plotter.draw()
            self.timeSlot = time.time()  # 更新时间戳
        else:
            # 绘制最终结果图
            if self.drawing != 0:
                if self.plotter is not None:  # 若绘制了动画，则保存并关闭动画
                    self.plotter.createAnimation()
                    self.plotter.close()
                if pop.ObjV.shape[1] == 2 or pop.ObjV.shape[1] == 3:
                    figureName = 'Pareto Front Plot'
                    self.plotter = ea.PointScatter(self.problem.M, grid=True, legend=True, title=figureName, saveName=self.dirName + figureName)
                    self.plotter.add(self.problem.ReferObjV, color='gray', alpha=0.1, label='True PF')
                    self.plotter.add(pop.ObjV, color='red', label='MOEA PF')
                    self.plotter.draw()
                else:
                    figureName = 'Parallel Coordinate Plot'
                    self.plotter = ea.ParCoordPlotter(self.problem.M, grid=True, legend=True, title=figureName, saveName=self.dirName + figureName)
                    self.plotter.add(self.problem.TinyReferObjV, color='gray', alpha=0.5, label='True Objective Value')
                    self.plotter.add(pop.ObjV, color='red', label='MOEA Objective Value')
                    self.plotter.draw()

    def stat(self, pop):

        """
        描述:
            该函数用于分析当代种群的信息。
            该函数会在terminated()函数里被调用。

        输入参数:
            pop : class <Population> - 种群对象。

        输出参数:
            无输出参数。

        """

        feasible = np.where(np.all(pop.CV <= 0, 1))[0] if pop.CV is not None else np.arange(pop.sizes)  # 找到满足约束条件的个体的下标
        if len(feasible) > 0:
            feasiblePop = pop[feasible]  # 获取满足约束条件的个体
            if self.logTras != 0 and self.currentGen % self.logTras == 0:
                self.logging(feasiblePop)  # 记录日志
                if self.verbose:
                    self.display()  # 打印日志
            self.draw(feasiblePop)  # 展示输出

    def terminated(self, pop):

        """
        描述:
            该函数用于判断是否应该终止进化，population为传入的种群对象。
            该函数会在各个具体的算法类的run()函数中被调用。

        输入参数:
            pop : class <Population> - 种群对象。

        输出参数:
            True / False。
            
        """

        self.check(pop)  # 检查种群对象的关键属性是否有误
        self.stat(pop)  # 进行统计分析
        self.passTime += time.time() - self.timeSlot  # 更新耗时
        # 调用outFunc()
        if self.outFunc is not None:
            if not callable(self.outFunc):
                raise RuntimeError('outFunc must be a function. (如果定义了outFunc，那么它必须是一个函数。)')
            self.outFunc(self, pop)
        self.timeSlot = time.time()  # 更新时间戳
        # 判断是否终止进化，
        if self.MAXGEN is None and self.MAXTIME is None and self.MAXEVALS is None:
            raise RuntimeError('error: MAXGEN, MAXTIME, and MAXEVALS cannot be all None. (MAXGEN, MAXTIME, 和MAXEVALS不能全为None)')
        terminatedFlag = False
        if self.MAXGEN is not None and self.currentGen + 1 >= self.MAXGEN:  # 由于代数是从0数起，因此在比较currentGen和MAXGEN时需要对currentGen加1
            self.stopMsg = 'The algotirhm stepped because it exceeded the generation limit.'
            terminatedFlag = True
        if self.MAXTIME is not None and self.passTime >= self.MAXTIME:
            self.stopMsg = 'The algotirhm stepped because it exceeded the time limit.'
            terminatedFlag = True
        if self.MAXEVALS is not None and self.evalsNum >= self.MAXEVALS:
            self.stopMsg = 'The algotirhm stepped because it exceeded the function evaluation limit.'
            terminatedFlag = True
        if terminatedFlag:
            return True
        else:
            self.currentGen += 1  # 进化代数+1
            return False

    def finishing(self, pop, globalNDSet=None):

        """
        描述:
            进化完成后调用的函数。

        输入参数:
            pop : class <Population> - 种群对象。
            
            globalNDSet : class <Population> - (可选参数)全局存档。

        输出参数:
            [NDSet, pop]，其中pop为种群类型；NDSet的类型与pop的一致。

        """
        if globalNDSet is None:
            # 得到非支配种群
            [levels, _] = ea.ndsortDED(pop.ObjV, needLevel=1, CV=pop.CV, maxormins=self.problem.maxormins)  # 非支配分层
            NDSet = pop[np.where(levels == 1)[0]]  # 只保留种群中的非支配个体，形成一个非支配种群
            if NDSet.CV is not None:  # CV不为None说明有设置约束条件
                NDSet = NDSet[np.where(np.all(NDSet.CV <= 0, 1))[0]]  # 最后要彻底排除非可行解
            # 排除重复解
            NDSet = NDSet[np.unique(NDSet.ObjV,return_index=True,axis=0)[1]]
        else:
            NDSet = globalNDSet
        if self.logTras != 0 and NDSet.sizes != 0 and (
                len(self.log['gen']) == 0 or self.log['gen'][-1] != self.currentGen):  # 补充记录日志和输出
            self.logging(NDSet)
            if self.verbose:
                self.display()
        self.passTime += time.time() - self.timeSlot  # 更新用时记录，因为已经要结束，因此不用再更新时间戳
        self.draw(NDSet, EndFlag=True)  # 显示最终结果图
        if self.plotter is not None:
            self.plotter.show()
        # 返回帕累托最优个体以及最后一代种群
        return [NDSet, pop]


class SoeaAlgorithm(Algorithm):  # 单目标优化算法类父类

    """
    描述:
        此为单目标进化优化算法类的父类，所有单目标优化算法类均继承自该父类。

    对比于父类该类新增的变量和函数:

        drawing         : int  - 绘图方式的参数，
                                 0表示不绘图；
                                 1表示绘制进化过程中种群的平均及最优目标函数值变化图；
                                 2表示实时绘制历史最优个体的目标空间过程动画；
                                 3表示实时绘制历史最优个体的决策空间动态图。

        trappedValue    : int  - 进化算法陷入停滞的判断阈值。

        maxTrappedCount : int  - “进化停滞”计数器最大上限值。
            
        ----------------- 以下为用户不需要设置的属性 -----------------

        BestIndi        : class <Population> - 存储算法所找到的最优的个体。

        trace           : dict - 进化记录器，可以看作是一个内部日志，用于记录每一代种群的一些信息。
                                 它与算法类的log类似，它有两个键：'f_best'以及'f_avg'。
                                 'f_best'的键值为一个list列表，存储着每一代种群最优个体的目标函数值；
                                 'f_avg'的键值为一个list列表，存储着每一代种群所有个体的平均目标函数值。

        trappedCount    : int  - “进化停滞”计数器。

        draw()          : 绘图函数。
        
    """

    def __init__(self,
                 problem,
                 population,
                 MAXGEN = None,
                 MAXTIME = None,
                 MAXEVALS = None,
                 MAXSIZE = None,
                 logTras = None,
                 verbose = None,
                 outFunc = None,
                 drawing = None,
                 trappedValue = None,
                 maxTrappedCount = None,
                 dirName = None,
                 **kwargs):

        """
        描述: 
            在该构造函数里只初始化静态参数以及对动态参数进行定义。
        
        """

        super().__init__(problem, population, MAXGEN, MAXTIME, MAXEVALS, MAXSIZE, logTras, verbose, outFunc, dirName)  # 先调用父类构造函数
        self.trappedValue = 0 if trappedValue is None else trappedValue
        self.maxTrappedCount = 1000 if maxTrappedCount is None else maxTrappedCount
        self.drawing = 1 if drawing is None else drawing
        # 以下为用户不需要设置的属性
        self.BestIndi = None  # 存储算法所找到的最优的个体
        self.trace = None  # 进化记录器
        self.trappedCount = None  # 定义trappedCount，在initialization()才对其进行初始化为0
        self.plotter = None  # 存储绘图对象

    def __str__(self):
        info = {}
        info['Algorithm Name'] = self.name
        info['Algorithm MAXGEN'] = self.MAXGEN
        info['Algorithm MAXTIME'] = self.MAXTIME
        info['Algorithm MAXEVALS'] = self.MAXEVALS
        info['Algorithm trappedValue'] = self.trappedValue
        info['Algorithm maxTrappedCount'] = self.maxTrappedCount
        return str(info)

    def initialization(self):

        """
        描述: 
            该函数用于在进化前对算法类的一些动态参数进行初始化操作。
            该函数需要在执行算法类的run()函数的一开始被调用，同时开始计时，
            以确保所有这些参数能够被正确初始化。
        
        """

        self.ax = None  # 初始化ax
        self.passTime = 0  # 初始化passTime
        self.trappedCount = 0  # 初始化“进化停滞”计数器
        self.currentGen = 0  # 初始为第0代
        self.evalsNum = 0  # 初始化评价次数为0
        self.BestIndi = ea.Population(None, None, 0)  # 初始化BestIndi为空的种群对象
        self.log = {'gen': [], 'eval': []} if self.logTras != 0 else None  # 初始化log
        self.trace = {'f_best': [], 'f_avg': []}  # 重置trace
        # 开始计时
        self.timeSlot = time.time()

    def logging(self, pop):

        """
        描述:
            用于在进化过程中记录日志。该函数在stat()函数里面被调用。
            如果需要在日志中记录其他数据，需要在自定义算法类中重写该函数。

        输入参数:
            pop : class <Population> - 种群对象。

        输出参数:
            无输出参数。

        """

        self.passTime += time.time() - self.timeSlot  # 更新用时记录，不计算logging的耗时
        if len(self.log['gen']) == 0:  # 初始化log的各个键值
            self.log['f_opt'] = []
            self.log['f_max'] = []
            self.log['f_avg'] = []
            self.log['f_min'] = []
            self.log['f_std'] = []
        self.log['gen'].append(self.currentGen)
        self.log['eval'].append(self.evalsNum)  # 记录评价次数
        self.log['f_opt'].append(self.BestIndi.ObjV[0][0])  # 记录算法所找到的最优个体的目标函数值
        self.log['f_max'].append(np.max(pop.ObjV))
        self.log['f_avg'].append(np.mean(pop.ObjV))
        self.log['f_min'].append(np.min(pop.ObjV))
        self.log['f_std'].append(np.std(pop.ObjV))
        self.timeSlot = time.time()  # 更新时间戳

    def draw(self, pop, EndFlag=False):

        """
        描述:
            该函数用于在进化过程中进行绘图。该函数在stat()以及finishing函数里面被调用。

        输入参数:
            pop     : class <Population> - 种群对象。
            
            EndFlag : bool - 表示是否是最后一次调用该函数。

        输出参数:
            无输出参数。

        """

        if not EndFlag:
            self.passTime += time.time() - self.timeSlot  # 更新用时记录，不计算画图的耗时
            # 绘制动画
            if self.drawing == 2:
                # 绘制迭代变化图
                trace_best = np.array(self.trace['f_best'])
                if len(trace_best) > 10:
                    xtickList = list(range(len(trace_best) - 10 + 1, len(trace_best) + 1))
                    trace_best = trace_best[len(trace_best) - 10:]
                else:
                    xtickList = list(range(1, len(trace_best)+1))
                num = len(trace_best)
                if self.plotter is None:
                    self.plotter = ea.ParCoordPlotter(num, grid=True, legend=True, title='Objective Value Trace Plot')
                self.plotter.refresh()
                self.plotter.Dimension = num
                self.plotter.xtickList = xtickList
                self.plotter.add(trace_best, color='blue', label='Best Objective Value')
                self.plotter.draw()
            elif self.drawing == 3:
                # 绘制决策空间动态图
                if self.plotter is None:
                    self.plotter = ea.ParCoordPlotter(self.problem.Dim, grid=True, legend=True, title='Variables Value Plot')
                self.plotter.refresh()
                self.plotter.add(pop.Phen, marker='o', color='blue', label='Variables Value at ' + str(self.currentGen) + ' Generation')
                self.plotter.draw()
            self.timeSlot = time.time()  # 更新时间戳
        else:
            # 绘制最终结果图
            if self.drawing != 0:
                if self.plotter is not None:  # 若绘制了动画，则保存并关闭动画
                    self.plotter.createAnimation()
                    self.plotter.close()
                trace_best = np.array(self.trace['f_best'])
                trace_avg = np.array(self.trace['f_avg'])
                num = len(trace_best)
                figureName = 'Trace Plot'
                self.plotter = ea.ParCoordPlotter(num, grid=True, legend=True, title=figureName, coordLabels=['Generation Number', 'Value'], saveName=self.dirName + figureName)
                self.plotter.add(trace_avg, color='blue', label='Average Objective Value')
                self.plotter.add(trace_best, color='orange', label='Best Objective Value')
                self.plotter.draw()

    def stat(self, pop):

        """
        描述:
            该函数用于分析、记录和打印当代种群的信息。
            该函数会在terminated()函数里被调用。

        输入参数:
            pop : class <Population> - 种群对象。

        输出参数:
            无输出参数。

        """

        # 进行进化记录
        feasible = np.where(np.all(pop.CV <= 0, 1))[0] if pop.CV is not None else np.arange(pop.sizes)  # 找到满足约束条件的个体的下标
        if len(feasible) > 0:
            feasiblePop = pop[feasible]
            bestIndi = feasiblePop[np.argmax(feasiblePop.FitnV)]  # 获取最优个体
            if self.BestIndi.sizes == 0:
                self.BestIndi = bestIndi  # 初始化global best individual
            else:
                delta = (self.BestIndi.ObjV - bestIndi.ObjV) * self.problem.maxormins if \
                    self.problem.maxormins is not None else self.BestIndi.ObjV - bestIndi.ObjV
                # 更新“进化停滞”计数器
                if np.abs(delta) < self.trappedValue:
                    self.trappedCount += 1
                else:
                    self.trappedCount = 0
                # 更新global best individual
                if delta > 0:
                    self.BestIndi = bestIndi
            # 更新trace
            self.trace['f_best'].append(bestIndi.ObjV[0][0])
            self.trace['f_avg'].append(np.mean(feasiblePop.ObjV))
            if self.logTras != 0 and self.currentGen % self.logTras == 0:
                self.logging(feasiblePop)  # 记录日志
                if self.verbose:
                    self.display()  # 打印日志
            self.draw(self.BestIndi)  # 展示输出

    def terminated(self, pop):

        """
        描述:
            该函数用于判断是否应该终止进化，population为传入的种群对象。
            该函数会在各个具体的算法类的run()函数中被调用。

        输入参数:
            pop : class <Population> - 种群对象。

        输出参数:
            True / False。

        """

        self.check(pop)  # 检查种群对象的关键属性是否有误
        self.stat(pop)  # 分析记录当代种群的数据
        self.passTime += time.time() - self.timeSlot  # 更新耗时
        # 调用outFunc()
        if self.outFunc is not None:
            if not callable(self.outFunc):
                raise RuntimeError('outFunc must be a function. (如果定义了outFunc，那么它必须是一个函数。)')
            self.outFunc(self, pop)
        self.timeSlot = time.time()  # 更新时间戳
        # 判断是否终止进化，
        if self.MAXGEN is None and self.MAXTIME is None and self.MAXEVALS is None:
            raise RuntimeError('error: MAXGEN, MAXTIME, and MAXEVALS cannot be all None. (MAXGEN, MAXTIME, 和MAXEVALS不能全为None)')
        terminatedFlag = False
        if self.MAXGEN is not None and self.currentGen + 1 >= self.MAXGEN:  # 由于代数是从0数起，因此在比较currentGen和MAXGEN时需要对currentGen加1
            self.stopMsg = 'The algotirhm stepped because it exceeded the generation limit.'
            terminatedFlag = True
        if self.MAXTIME is not None and self.passTime >= self.MAXTIME:
            self.stopMsg = 'The algotirhm stepped because it exceeded the time limit.'
            terminatedFlag = True
        if self.MAXEVALS is not None and self.evalsNum >= self.MAXEVALS:
            self.stopMsg = 'The algotirhm stepped because it exceeded the function evaluation limit.'
            terminatedFlag = True
        if self.trappedCount >= self.maxTrappedCount:
            self.stopMsg = 'The algotirhm stepped because it exceeded the trapped count limit.'
            terminatedFlag = True
        if terminatedFlag:
            return True
        else:
            self.currentGen += 1  # 进化代数+1
            return False

    def finishing(self, pop):

        """
        描述:
            进化完成后调用的函数。

        输入参数:
            pop : class <Population> - 种群对象。

        输出参数:
            [self.BestIndi, pop]，其中pop为种群类型；BestIndi的类型与pop的一致。
        
        注意:
            若没有找到可行解，则返回的self.BestIndi为None。
        
        """

        feasible = np.where(np.all(pop.CV <= 0, 1))[0] if pop.CV is not None else np.arange(pop.sizes)  # 找到满足约束条件的个体的下标
        if len(feasible) > 0:
            feasiblePop = pop[feasible]
            if self.logTras != 0 and (len(self.log['gen']) == 0 or self.log['gen'][-1] != self.currentGen):  # 补充记录日志和输出
                self.logging(feasiblePop)
                if self.verbose:
                    self.display()
        self.passTime += time.time() - self.timeSlot  # 更新用时记录，因为已经要结束，因此不用再更新时间戳
        self.draw(pop, EndFlag=True)  # 显示最终结果图
        if self.plotter:
            self.plotter.show()
        # 返回最优个体以及最后一代种群
        return [self.BestIndi, pop]

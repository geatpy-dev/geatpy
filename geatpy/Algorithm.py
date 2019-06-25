# -*- coding: utf-8 -*-

class Algorithm:
    
    """
Algorithm : class - 算法模板类，所有算法模板都要继承该父类。

描述:
    算法设置类是用来存储与算法运行参数设置相关信息的一个类。

属性:
    name          : str      - 算法名称（可以自由设置名称）。
    
    problem       : class <Problem> - 问题类的对象。
    
    MAXGEN        : int      - 最大进化代数。
    
    currentGen    : int      - 当前进化的代数。
    
    MAXTIME       : float    - 时间限制（单位：秒）。
    
    timeSlot      : flot     - 时间戳（单位：秒）。
    
    passTime      : float    - 已用时间（单位：秒）。
    
    MAXEVALS      : int      - 最大评价次数。
    
    evalsNum      : int      - 当前评价次数。
    
    MAXSIZE       : int      - 最优解的最大数目。
    
    population    : class <Population> - 种群对象。
    
    selFunc       : str      - 选择算子的名称。
    
    recFunc       : str      - 重组算子的名称。
    
    mutFunc       : str      - 变异算子的名称。
    
    drawing       : int      - 绘图方式的参数，0表示不绘图，1表示绘图，2表示实时绘制动态图。

函数:
    terminated()  : 根据参数的设置，计算是否需要终止进化，需要在继承类中实现。
    
    run() : 执行函数，需要在继承类中实现。
    
"""

    def __init__(self):
        self.name = None
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
        self.selFunc = None
        self.recFunc = None
        self.mutFunc = None
        self.drawing = None
    
    def iteration(self):
        pass
    
    def run(self):
        pass
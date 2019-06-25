# -*- coding: utf-8 -*-

class Problem:
    
    """
Problem : Class - 问题类

描述:
    问题类是用来存储与待求解问题相关信息的一个类。

属性:
    name      : str   - 问题名称（可以自由设置名称）。
    
    M         : int   - 目标维数，即有多少个优化目标。
    
    maxormins : list  - 目标最小最大化标记列表，1表示最小化，-1表示最大化，例如：
                        [1,1,-1,-1]，表示前2个目标是最小化，后2个目标是最大化。
                           
    Dim       : int   - 决策变量维数，即有多少个决策变量。
    
    varTypes  : array - 连续或离散标记行向量，0表示对应的决策变量是连续的；
                        1表示对应的变量是离散的。
    
    ranges    : array - 决策变量范围矩阵，第一行对应决策变量的下界，第二行对应决策变量的上界。
    
    borders   : array - 决策变量范围的边界矩阵，第一行对应决策变量的下边界，第二行对应决策变量的上边界，
                        0表示范围中不含边界，1表示范围包含边界。

函数:
    aimFuc() : 目标函数，需要在继承类中实现。
    
    calBest() : 计算理论最优值的函数，需要在继承类中实现。

"""

    def __init__(self):
        self.name = None
        self.M = None
        self.maxormins = None
        self.Dim = None
        self.varTypes = None
        self.ranges = None
        self.borders = None
    
    def aimFuc(self):
        pass
    
    def calBest(self):
        pass
    

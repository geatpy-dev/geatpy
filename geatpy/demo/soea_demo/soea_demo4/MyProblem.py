# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

"""
该案例展示了一个利用单目标进化算法实现句子匹配的应用实例。
"""

class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        # 定义需要匹配的句子
        strs = 'Tom is a little boy, isn\'t he? Yes he is, he is a good and smart child and he is always ready to help others, all in all we all like him very much.'
        self.words = []
        for c in strs:
            self.words.append(ord(c)) # 把字符串转成ASCII码
        M = 1 # 初始化M（目标维数）
        maxormins = [1] # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = len(self.words) # 初始化Dim（决策变量维数）
        varTypes = [1] * Dim # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [32] * Dim # 决策变量下界
        ub = [122] * Dim # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界
        ubin = [1] * Dim # 决策变量上边界
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
        
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        diff = np.sum((Vars - self.words)**2, 1)
        pop.ObjV = np.array([diff]).T # 把求得的目标函数值赋值给种群pop的ObjV
    
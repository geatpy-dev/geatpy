# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

"""
该案例展示了一个利用单目标进化算法实现句子匹配的应用实例。
"""

class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self):
        self.name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        # 定义需要匹配的句子
        strs = 'Tom is a little boy, isn\'t he? Yes he is, he is a good and smart child and he is always ready to help others, all in all we all like him very much.'
        self.words = []
        for c in strs:
            self.words.append(ord(c)) # 把字符串转成ASCII码
        self.M = 1 # 初始化M（目标维数）
        self.maxormins = [1] # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = len(self.words) # 初始化Dim（决策变量维数）
        self.varTypes = np.array([1] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [32] * self.Dim # 决策变量下界
        ub = [122] * self.Dim # 决策变量上界
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
        
    def aimFuc(self, Vars, CV): # 目标函数设置
        diff = np.sum((Vars - self.words)**2, 1)
        return [np.array([diff]).T, CV]
    
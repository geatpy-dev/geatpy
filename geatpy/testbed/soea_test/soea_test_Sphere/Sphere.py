# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class Sphere(ea.Problem): # 继承Problem父类
    def __init__(self, D = 2):
        self.name = 'Beale' # 初始化name（函数名称，可以随意设置）
        self.M = 1 # 初始化M（目标维数）
        self.maxormins = [1] * self.M # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = D # 初始化Dim（决策变量维数）
        self.varTypes = np.array([0] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [-5.12] * self.Dim # 决策变量下界
        ub = [5.12] * self.Dim # 决策变量上界
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
    
    def aimFuc(self, Vars, CV):
        f = np.array([np.sum(Vars**2, 1)]).T
        return f, CV
    
    def calBest(self):
        realBestObjV = np.array([[0]])
        return realBestObjV
    
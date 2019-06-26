# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class Shubert(ea.Problem): # 继承Problem父类
    def __init__(self, M = 1, D = 2):
        self.name = 'Shubert' # 初始化name（函数名称，可以随意设置）
        self.M = M # 初始化M（目标维数）
        self.maxormins = [1] * self.M # 初始化maxormins（目标最小最大化标记列表）
        self.Dim = D # 初始化Dim（决策变量维数）
        self.varTypes = np.array([0] * self.Dim) # 初始化varTypes（决策变量的类型）
        lb = [-10] * self.Dim # 决策变量下界
        ub = [10] * self.Dim # 决策变量上界
        self.ranges = np.array([lb, ub]) # 初始化ranges（决策变量范围矩阵）
        lbin = [1] * self.Dim
        ubin = [1] * self.Dim
        self.borders = np.array([lbin, ubin]) # 初始化borders（决策变量范围边界矩阵）
    
    def aimFuc(self, Vars, CV):
        x = Vars[:, [0]]
        y = Vars[:, [1]]
        f = ((1*np.cos((1+1)*x+1)) + (2*np.cos((2+1)*x+2)) + (3*np.cos((3+1)*x+3)) + 
             (4*np.cos((4+1)*x+4)) + (5*np.cos((5+1)*x+5))) * ((1*np.cos((1+1)*y+1)) + 
             (2*np.cos((2+1)*y+2)) + (3*np.cos((3+1)*y+3)) + (4*np.cos((4+1)*y+4)) + 
             (5*np.cos((5+1)*y+5)))
        return f, CV
    
    def calBest(self):
        realBestObjV = np.array([[-186.731]])
        
        return realBestObjV
    
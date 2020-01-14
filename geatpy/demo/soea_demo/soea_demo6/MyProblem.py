# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea
from sklearn import svm
from sklearn.model_selection import cross_val_score
from multiprocessing.dummy import Pool as ThreadPool

"""
该案例展示了如何利用进化算法优化SVM中的两个参数：C和Gamma。
在执行本案例前，需要确保正确安装sklearn，以保证SVM部分的代码能够正常执行。
本函数需要用到一个外部数据集，存放在同目录下的iris.data中。
有关该数据集的详细描述详见http://archive.ics.uci.edu/ml/datasets/Iris
"""

class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        M = 1 # 初始化M（目标维数）
        maxormins = [-1] # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 2 # 初始化Dim（决策变量维数）
        varTypes = [0, 0] # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [2**(-8)] * Dim # 决策变量下界
        ub = [2**8] * Dim # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1] * Dim # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
        # 目标函数计算中用到的一些数据
        fp = open('iris.data')
        datas = []
        data_targets = []
        for line in fp.readlines():
            line_data = line.strip('\n').split(',')
            data = []
            for i in line_data[0:4]:
                data.append(float(i))
            data_targets.append(line_data[4])
            datas.append(data)
        self.dataTarget = np.array(data_targets)
        self.data = np.array(datas)
        fp.close()
    
    def aimFunc(self, pop): # 目标函数，采用多线程加速计算
        Vars = pop.Phen # 得到决策变量矩阵
        pop.ObjV = np.zeros((pop.sizes, 1)) # 初始化种群个体目标函数值列向量
        def subAimFunc(i):
            C = Vars[i, 0]
            G = Vars[i, 1]
            svc = svm.SVC(C=C, kernel='rbf', gamma=G).fit(self.data, self.dataTarget)
            scores = cross_val_score(svc, self.data, self.dataTarget, cv=10) # 计算交叉验证的得分
            pop.ObjV[i] = scores.mean() # 把交叉验证的平均得分作为目标函数值
        pool = ThreadPool(2) # 设置池的大小
        pool.map(subAimFunc, list(range(pop.sizes)))
    
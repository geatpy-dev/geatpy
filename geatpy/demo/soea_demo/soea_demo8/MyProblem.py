# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import geatpy as ea


class MyProblem(ea.Problem):  # 继承Problem父类
    def __init__(self):
        # 目标函数计算中用到的一些数据
        self.datas = np.loadtxt('data.csv', delimiter=',')  # 读取数据
        self.k = 4  # 分类数目
        # 问题类设置
        name = 'MyProblem'  # 初始化name（函数名称，可以随意设置）
        M = 1  # 初始化M（目标维数）
        maxormins = [1]  # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = self.datas.shape[1] * self.k  # 初始化Dim
        varTypes = [0] * Dim  # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = list(np.min(self.datas, 0)) * self.k  # 决策变量下界
        ub = list(np.max(self.datas, 0)) * self.k  # 决策变量上界
        lbin = [1] * Dim  # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1] * Dim  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)

    def aimFunc(self, pop):  # 目标函数
        centers = pop.Phen.reshape(int(pop.sizes * self.k), int(pop.Phen.shape[1] / self.k))  # 得到聚类中心
        dis = ea.cdist(centers, self.datas, 'euclidean')  # 计算距离
        dis_split = dis.reshape(pop.sizes, self.k, self.datas.shape[0])  # 分割距离矩阵，把各个聚类中心到各个点之间的距离的数据分开
        labels = np.argmin(dis_split, 1)[0]  # 得到聚类标签值
        uni_labels = np.unique(labels)
        for i in range(len(uni_labels)):
            centers[uni_labels[i], :] = np.mean(self.datas[np.where(labels == uni_labels[i])[0], :], 0)
        # 直接修改染色体为已知的更优值，加快收敛
        pop.Chrom = centers.reshape(pop.sizes, self.k * centers.shape[1])
        pop.Phen = pop.decoding()  # 染色体解码（要同步修改Phen，否则后面会导致数据不一致）
        dis = ea.cdist(centers, self.datas, 'euclidean')
        dis_split = dis.reshape(pop.sizes, self.k, self.datas.shape[0])
        pop.ObjV = np.sum(np.min(dis_split, 1), 1, keepdims=True)  # 计算个体的目标函数值

    def draw(self, centers):  # 绘制聚类效果图
        dis = ea.cdist(centers, self.datas, 'euclidean')
        dis_split = dis.reshape(1, self.k, self.datas.shape[0])
        labels = np.argmin(dis_split, 1)[0]
        colors = ['r', 'g', 'b', 'y']
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for i in range(self.k):
            idx = np.where(labels == i)[0]  # 找到同一类的点的下标
            datas = self.datas[idx, :]
            ax.scatter(datas[:, 0], datas[:, 1], datas[:, 2], c=colors[i])

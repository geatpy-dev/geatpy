# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea
from sklearn import svm
from sklearn import preprocessing
from sklearn.model_selection import cross_val_score
import multiprocessing as mp
from multiprocessing import Pool as ProcessPool
from multiprocessing.dummy import Pool as ThreadPool

"""
该案例展示了如何利用进化算法+多进程/多线程来优化SVM中的两个参数：C和Gamma。
在执行本案例前，需要确保正确安装sklearn，以保证SVM部分的代码能够正常执行。
本函数需要用到一个外部数据集，存放在同目录下的iris.data中，
并且把iris.data按3:2划分为训练集数据iris_train.data和测试集数据iris_test.data。
有关该数据集的详细描述详见http://archive.ics.uci.edu/ml/datasets/Iris
在执行脚本main.py中设置PoolType字符串来控制采用的是多进程还是多线程。
注意：使用多进程时，程序必须以“if __name__ == '__main__':”作为入口，
      这个是multiprocessing的多进程模块的硬性要求。
"""


class MyProblem(ea.Problem):  # 继承Problem父类
    def __init__(self, PoolType):  # PoolType是取值为'Process'或'Thread'的字符串
        name = 'MyProblem'  # 初始化name（函数名称，可以随意设置）
        M = 1  # 初始化M（目标维数）
        maxormins = [-1]  # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 2  # 初始化Dim（决策变量维数）
        varTypes = [0, 0]  # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [2 ** (-8)] * Dim  # 决策变量下界
        ub = [2 ** 8] * Dim  # 决策变量上界
        lbin = [1] * Dim  # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1] * Dim  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
        # 目标函数计算中用到的一些数据
        fp = open('iris_train.data')
        datas = []
        data_targets = []
        for line in fp.readlines():
            line_data = line.strip('\n').split(',')
            data = []
            for i in line_data[0:4]:
                data.append(float(i))
            datas.append(data)
            data_targets.append(line_data[4])
        fp.close()
        self.data = preprocessing.scale(np.array(datas))  # 训练集的特征数据（归一化）
        self.dataTarget = np.array(data_targets)
        # 设置用多线程还是多进程
        self.PoolType = PoolType
        if self.PoolType == 'Thread':
            self.pool = ThreadPool(2)  # 设置池的大小
        elif self.PoolType == 'Process':
            num_cores = int(mp.cpu_count())  # 获得计算机的核心数
            self.pool = ProcessPool(num_cores)  # 设置池的大小

    def aimFunc(self, pop):  # 目标函数，采用多线程加速计算
        Vars = pop.Phen  # 得到决策变量矩阵
        args = list(
            zip(list(range(pop.sizes)), [Vars] * pop.sizes, [self.data] * pop.sizes, [self.dataTarget] * pop.sizes))
        if self.PoolType == 'Thread':
            pop.ObjV = np.array(list(self.pool.map(subAimFunc, args)))
        elif self.PoolType == 'Process':
            result = self.pool.map_async(subAimFunc, args)
            result.wait()
            pop.ObjV = np.array(result.get())

    def test(self, C, G):  # 代入优化后的C、Gamma对测试集进行检验
        # 读取测试集数据
        fp = open('iris_test.data')
        datas = []
        data_targets = []
        for line in fp.readlines():
            line_data = line.strip('\n').split(',')
            data = []
            for i in line_data[0:4]:
                data.append(float(i))
            datas.append(data)
            data_targets.append(line_data[4])
        fp.close()
        data_test = preprocessing.scale(np.array(datas))  # 测试集的特征数据（归一化）
        dataTarget_test = np.array(data_targets)  # 测试集的标签数据
        svc = svm.SVC(C=C, kernel='rbf', gamma=G).fit(self.data, self.dataTarget)  # 创建分类器对象并用训练集的数据拟合分类器模型
        dataTarget_predict = svc.predict(data_test)  # 采用训练好的分类器对象对测试集数据进行预测
        print("测试集数据分类正确率 = %s%%" % (
                len(np.where(dataTarget_predict == dataTarget_test)[0]) / len(dataTarget_test) * 100))


def subAimFunc(args):
    i = args[0]
    Vars = args[1]
    data = args[2]
    dataTarget = args[3]
    C = Vars[i, 0]
    G = Vars[i, 1]
    svc = svm.SVC(C=C, kernel='rbf', gamma=G).fit(data, dataTarget)  # 创建分类器对象并用训练集的数据拟合分类器模型
    scores = cross_val_score(svc, data, dataTarget, cv=30)  # 计算交叉验证的得分
    ObjV_i = [scores.mean()]  # 把交叉验证的平均得分作为目标函数值
    return ObjV_i

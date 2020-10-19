# -*- coding: utf-8 -*-
import numpy as np
import xlrd
import geatpy as ea
from sklearn import svm
from sklearn import preprocessing
from sklearn.model_selection import cross_val_score
from scoop import futures

"""
该案例展示了如何利用SCOOP库进行分布式加速计算Geatpy进化算法程序，
本案例和soea_demo6类似，同样是用进化算法来优化SVM的参数C和Gamma，
不同的是，本案例选用更庞大的数据集，使得每次训练SVM模型时耗时更高，从而更适合采用分布式加速计算。
该数据集存放在同目录下的Data_User_Modeling_Dataset_Hamdi Tolga KAHRAMAN.xls中，
有关该数据集的详细描述详见http://archive.ics.uci.edu/ml/datasets/User+Knowledge+Modeling。
在执行本案例前，需要确保正确安装sklearn以及SCOOP，以保证SVM和SCOOP部分的代码能够正常执行。
SCOOP安装方法：控制台执行命令pip install scoop
分布式加速计算注意事项：
1.当aimFunc()函数十分耗时，比如无法矩阵化计算、或者是计算单个个体的目标函数值就需要很长时间时，
  适合采用分布式计算，否则贸然采用分布式计算反而会大大降低性能。
2.分布式执行方法：python -m scoop -n 10 main.py 其中10表示把计算任务分发给10个workers。
  非分布式执行方法：python main.py
"""


class MyProblem(ea.Problem):  # 继承Problem父类
    def __init__(self):
        name = 'MyProblem'  # 初始化name（函数名称，可以随意设置）
        M = 1  # 初始化M（目标维数）
        maxormins = [-1]  # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 2  # 初始化Dim（决策变量维数）
        varTypes = [0, 0]  # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [2 ** (-8), 2 ** (-8)]  # 决策变量下界
        ub = [2 ** 8, 1]  # 决策变量上界
        lbin = [1] * Dim  # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1] * Dim  # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
        # 目标函数计算中用到的一些数据
        workbook = xlrd.open_workbook(
            "Data_User_Modeling_Dataset_Hamdi Tolga KAHRAMAN.xls")  # 打开文件，获取excel文件的workbook（工作簿）对象
        worksheet = workbook.sheet_by_name("Training_Data")  # 通过sheet名获得sheet对象
        self.data = np.vstack([worksheet.col_values(0)[1:],
                               worksheet.col_values(1)[1:],
                               worksheet.col_values(2)[1:],
                               worksheet.col_values(3)[1:],
                               worksheet.col_values(4)[1:]]).T  # 获取特征数据
        self.data = preprocessing.scale(self.data)  # 归一化特征数据
        self.dataTarget = worksheet.col_values(5)[1:]  # 获取标签数据

    def aimFunc(self, pop):  # 目标函数
        Vars = pop.Phen  # 得到决策变量矩阵
        args = list(
            zip(list(range(pop.sizes)), [Vars] * pop.sizes, [self.data] * pop.sizes, [self.dataTarget] * pop.sizes))
        pop.ObjV = np.array(list(futures.map(subAimFunc, args)))  # 调用SCOOP的map函数进行分布式计算，并构造种群所有个体的目标函数值矩阵ObjV

    def test(self, C, G):  # 代入优化后的C、Gamma对测试集进行检验
        # 读取测试集数据
        workbook = xlrd.open_workbook(
            "Data_User_Modeling_Dataset_Hamdi Tolga KAHRAMAN.xls")  # 打开文件，获取excel文件的workbook（工作簿）对象
        worksheet = workbook.sheet_by_name("Test_Data")  # 通过sheet名获得sheet对象
        data_test = np.vstack([worksheet.col_values(0)[1:],
                               worksheet.col_values(1)[1:],
                               worksheet.col_values(2)[1:],
                               worksheet.col_values(3)[1:],
                               worksheet.col_values(4)[1:]]).T  # 获取特征数据
        data_test = preprocessing.scale(data_test)  # 归一化特征数据
        dataTarget_test = worksheet.col_values(5)[1:]  # 获取标签数据
        svc = svm.SVC(C=C, kernel='rbf', gamma=G).fit(self.data, self.dataTarget)  # 创建分类器对象并用训练集的数据拟合分类器模型
        dataTarget_predict = svc.predict(data_test)  # 采用训练好的分类器对象对测试集数据进行预测
        print("测试集数据分类正确率 = %s%%" % (
                len(np.where(dataTarget_predict == dataTarget_test)[0]) / len(dataTarget_test) * 100))


def subAimFunc(args):  # 单独计算单个个体的目标函数值
    i = args[0]
    Vars = args[1]
    data = args[2]
    dataTarget = args[3]
    C = Vars[i, 0]
    G = Vars[i, 1]
    svc = svm.SVC(C=C, kernel='rbf', gamma=G).fit(data, dataTarget)  # 创建分类器对象并用训练集的数据拟合分类器模型
    scores = cross_val_score(svc, data, dataTarget, cv=20)  # 计算交叉验证的得分
    ObjV_i = [scores.mean()]  # 把交叉验证的平均得分作为目标函数值
    return ObjV_i

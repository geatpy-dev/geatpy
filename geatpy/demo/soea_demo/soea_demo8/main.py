# -*- coding: utf-8 -*-
from MyProblem import MyProblem  # 导入自定义问题接口
import geatpy as ea  # import geatpy

"""
该案例展示了如何利用进化算法进行仿k-means聚类（可称之为EA-KMeans算法）。问题的定义详见MyProblem.py。
本案例采用与k-means类似的聚类方法，采用展开的聚类中心点坐标作为染色体的编码，基本流程大致如下：
1) 初始化种群染色体。
2) 迭代进化（循环第3步至第6步），直到满足终止条件。
3) 重组变异，然后根据得到的新染色体计算出对应的聚类中心点。
4) 计算各数据点到聚类中心点的欧式距离。
5) 把与各中心点关联的数据点的坐标平均值作为新的中心点，并以此更新种群的染色体。
6) 把各中心点到与其关联的数据点之间的距离之和作为待优化的目标函数值。
注意：导入的数据是以列为特征的，即每一列代表一个特征（如第一列代表x，第二列代表y......）。
"""

if __name__ == '__main__':
    # 实例化问题对象
    problem = MyProblem()
    # 构建算法
    algorithm = ea.soea_DE_rand_1_bin_templet(problem,
                                              ea.Population(Encoding='RI', NIND=2),
                                              MAXGEN=20,  # 最大进化代数。
                                              logTras=1,  # 表示每隔多少代记录一次日志信息，0表示不记录。
                                              trappedValue=1e-4,  # 单目标优化陷入停滞的判断阈值。
                                              maxTrappedCount=10)  # 进化停滞计数器最大上限值。
    # 求解
    res = ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=True, drawLog=False, saveFlag=True)
    # 检验结果
    if res['success']:
        print('最优的聚类中心为：')
        Vars = res['Vars'][0, :]
        centers = Vars.reshape(problem.k, int(len(Vars) / problem.k))  # 得到最优的聚类中心
        print(centers)
        """=================================检验结果==============================="""
        problem.draw(centers)

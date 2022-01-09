# -*- coding: utf-8 -*-
from MyProblem import MyProblem  # 导入自定义问题接口
import geatpy as ea  # import geatpy

"""
本案例和soea_demo6类似，同样是用进化算法来搜索SVM的参数C和Gamma的最佳值。不同的是这里采用SCOOP来
在执行本案例前，需要确保正确安装sklearn以及SCOOP，以保证SVM和SCOOP部分的代码能够正常执行。
SCOOP安装方法：控制台执行命令pip install scoop
分布式加速计算注意事项：
1.当目标函数的计算十分耗时，比如计算单个个体的目标函数值就需要很长时间时，
  适合采用分布式计算，否则贸然采用分布式计算反而会大大降低性能。
2.分布式执行方法：python -m scoop -n 10 main.py 其中10表示把计算任务分发给10个workers。
  非分布式执行方法：python main.py
"""

if __name__ == '__main__':
    # 实例化问题对象
    problem = MyProblem()
    # 构建算法
    algorithm = ea.soea_DE_rand_1_bin_templet(problem,
                                              ea.Population(Encoding='RI', NIND=20),
                                              MAXGEN=30,  # 最大进化代数。
                                              logTras=1,  # 表示每隔多少代记录一次日志信息，0表示不记录。
                                              trappedValue=1e-6,  # 单目标优化陷入停滞的判断阈值。
                                              maxTrappedCount=10)  # 进化停滞计数器最大上限值。
    # 求解
    res = ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=True, drawLog=False, saveFlag=True)
    # 检验结果
    problem.test(C=res['Vars'][0, 0], G=res['Vars'][0, 1])

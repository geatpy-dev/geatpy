"""最小化目标函数值问题求解程序执行脚本main.py"""
import numpy as np
import geatpy as gea # 导入geatpy库
from mintemp1 import mintemp1 # 导入自定义的编程模板
from aimfuc import aimfuc # 导入自定义的目标函数接口
from punishing import punishing # 导入自定义的罚函数接口

# 获取函数接口地址
AIM_M = __import__('aimfuc')
PUN_M = __import__('punishing')
"""============================变量设置============================"""
x1 = [-5, 5]; x2 = [2, 10]          # 自变量的范围
b1 = [1, 1] ; b2 = [0, 1]           # 自变量的边界
ranges=np.vstack([x1, x2]).T        # 生成自变量的范围矩阵
borders=np.vstack([b1, b2]).T       # 生成自变量的边界矩阵
"""========================遗传算法参数设置========================="""
NIND = 10;               # 种群规模
MAXGEN = 50;             # 最大遗传代数
GGAP = 0.8;              # 代沟：子代与父代的重复率为(1-GGAP)
selectStyle = 'rws';     # 遗传算法的选择方式设为"rws"——轮盘赌选择
recombinStyle = 'xovdp'  # 遗传算法的重组方式，设为两点交叉
recopt = 0.9;            # 交叉概率
pm = 0.01;               # 变异概率
SUBPOP = 1               # 设置种群数为1
maxormin = 1             # 设置标记表明这是最小化目标
"""=======================调用编程模板进行种群进化==================="""
# 调用编程模板进行种群进化，得到种群进化和变量的追踪器以及运行时间
[pop_trace, var_trace, times] = mintemp1(AIM_M, 'aimfuc', PUN_M, 'punishing', ranges, borders, MAXGEN, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, maxormin)
"""=========================绘图及输出结果========================="""
# 传入pop_trace进行绘图
gea.trcplot(pop_trace, [['种群最优个体目标函数值'], ['种群个体平均适应度值', '种群最优个体适应度值']])
# 输出结果
best_gen = np.argmin(pop_trace[:, 0]) # 记录最优种群是在哪一代
print('最优的目标函数值为：', np.min(pop_trace[:, 0]))
print('最优的控制变量值为：')
for i in range(var_trace.shape[1]):
	print(var_trace[best_gen, i])
print('最优的一代是第',best_gen + 1,'代')
print('用时：', times, '秒')
"""
执行脚本main.py
描述：
	该demo是展示如何计算带约束的单目标优化问题
	本案例通过自定义算法模板"mintemp1.py"来解决该问题
	其中目标函数和约束条件写在aimfuc.py文件中，适应度罚函数写在罚函数文件punishing.py中
"""
import numpy as np
import geatpy as gea # 导入geatpy库
from mintemp1 import mintemp1 # 导入自定义的编程模板

# 获取函数接口地址
AIM_M = __import__('aimfuc')        # 目标函数
PUN_M = __import__('punishing')     # 罚函数
"""============================变量设置============================"""
x1 = [-5, 5]; x2 = [2, 10]          # 自变量的范围
b1 = [1, 1] ; b2 = [0, 1]           # 自变量的边界
ranges=np.vstack([x1, x2]).T        # 生成自变量的范围矩阵
borders=np.vstack([b1, b2]).T       # 生成自变量的边界矩阵
"""========================遗传算法参数设置========================="""
NIND = 50                # 种群规模
MAXGEN = 500              # 最大遗传代数
GGAP = 0.8               # 代沟：子代与父代的重复率为(1-GGAP)
selectStyle = 'rws'      # 遗传算法的选择方式设为"rws"——轮盘赌选择
recombinStyle = 'xovdp'  # 遗传算法的重组方式，设为两点交叉
recopt = 0.9             # 交叉概率
pm = 0.01                # 变异概率
SUBPOP = 1               # 设置种群数为1
maxormin = 1             # 设置标记表明这是最小化目标
"""=======================调用编程模板进行种群进化==================="""
# 调用编程模板进行种群进化，得到种群进化和变量的追踪器以及运行时间
[pop_trace, var_trace, times] = mintemp1(AIM_M, 'aimfuc', PUN_M, 'punishing', ranges, borders, MAXGEN, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, maxormin)
"""=========================绘图及输出结果========================="""
# 传入pop_trace进行绘图
gea.trcplot(pop_trace, [['各代种群最优目标函数值'], ['各代种群个体平均适应度值', '各代种群最优个体适应度值']], ['demo_result1', 'demo_result2'])
# 输出结果
best_gen = np.argmin(pop_trace[:, 0]) # 记录最优种群是在哪一代
print('最优的目标函数值为：', np.min(pop_trace[:, 0]))
print('最优的控制变量值为：')
for i in range(var_trace.shape[1]):
	print(var_trace[best_gen, i])
print('最优的一代是第',best_gen + 1,'代')
print('用时：', times, '秒')
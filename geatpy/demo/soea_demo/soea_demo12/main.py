# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea

class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'Shortest_Path' # 初始化name（函数名称，可以随意设置）
        M = 1 # 初始化M（目标维数）
        maxormins = [1] # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 10 # 初始化Dim（决策变量维数）
        varTypes = [1] * Dim # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [0] * Dim # 决策变量下界
        ub = [9] * Dim # 决策变量上界
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub)
        # 设置每一个结点下一步可达的结点（结点从1开始数，因此列表nodes的第0号元素设为空列表表示无意义）
        self.nodes = [[], [2,3], [3,4,5], [5,6], [7,8], [4,6], [7,9], [8,9], [9,10], [10]]
        # 设置有向图中各条边的权重
        self.weights = {'(1, 2)':36, '(1, 3)':27, '(2, 4)':18, '(2, 5)':20, '(2, 3)':13, '(3, 5)':12, '(3, 6)':23,
                         '(4, 7)':11, '(4, 8)':32, '(5, 4)':16, '(5, 6)':30, '(6, 7)':12, '(6, 9)':38, '(7, 8)':20,
                         '(7, 9)':32, '(8, 9)':15, '(8, 10)':24, '(9, 10)':13}
    
    def decode(self, priority): # 将优先级编码的染色体解码得到一条从节点1到节点10的可行路径
        edges = [] # 存储边
        path = [1] # 结点1是路径起点
        while not path[-1] == 10: # 开始从起点走到终点
            currentNode = path[-1] # 得到当前所在的结点编号
            nextNodes = self.nodes[currentNode] # 获取下一步可达的结点组成的列表
            chooseNode = nextNodes[np.argmax(priority[np.array(nextNodes) - 1])] # 从NextNodes中选择优先级更高的结点作为下一步要访问的结点，因为结点从1数起，而下标从0数起，因此要减去1
            path.append(chooseNode)
            edges.append((currentNode, chooseNode))
        return path, edges

    def evalVars(self, Vars): # 目标函数
        N = Vars.shape[0]
        ObjV = np.zeros((N, 1)) # 初始化ObjV
        for i in range(N): # 遍历种群的每个个体，分别计算各个个体的目标函数值
            priority = Vars[i, :]
            path, edges = self.decode(priority) # 将优先级编码的染色体解码得到访问路径及经过的边
            pathLen = 0
            for edge in edges:
                key = str(edge) # 根据路径得到键值，以便根据键值找到路径对应的长度
                if not key in self.weights:
                    raise RuntimeError("Error in aimFunc: The path is invalid. (当前路径是无效的。)", path)
                pathLen += self.weights[key] # 将该段路径长度加入
            ObjV[i] = pathLen # 计算目标函数值，赋值给pop种群对象的ObjV属性
        return ObjV
## 执行脚本
if __name__ == "__main__":
    # 实例化问题对象
    problem = MyProblem()
    # 构建算法
    algorithm = ea.soea_EGA_templet(problem,
                                    ea.Population(Encoding='RI', NIND=4),
                                    MAXGEN=10,  # 最大进化代数
                                    logTras=1)  # 表示每隔多少代记录一次日志信息
    # 求解
    res = ea.optimize(algorithm, verbose=True, drawing=1, outputMsg=False, drawLog=False, saveFlag=True, dirName='result')
    print('最短路程为：%s'%(res['ObjV'][0][0]))
    print('最佳路线为：')
    best_journey, edges = problem.decode(res['Vars'][0])
    for i in range(len(best_journey)):
        print(int(best_journey[i]), end = ' ')
    print()

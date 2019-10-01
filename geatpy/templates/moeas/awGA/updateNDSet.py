# -*- coding: utf-8 -*-

import numpy as np
import geatpy as ea

def updateNDSet(population, maxormins, MAXSIZE, NDSet = None):
    
    """
描述:
    用于计算种群个体的适应度及更新全局非支配种群。
输入参数:
    population : Population - 种群对象。
    
    maxormins  : list - 优化目标的最大最小化标记列表。
    
    MAXSIZE    : int - 示全局非支配个体的最大数目。
    
    NDSet      : Population - (可选参数)全局非支配个体，
                              若缺省或为None时，NDSet为所传入种群的非支配个体组成的种群。

输出参数:
    NDSet      : Population - (可选参数)全局非支配种群。
    
    种群适应度FitnV已经在函数中进行更新，因此这该参数不用返回。
    """
    
    ObjV = maxormins * population.ObjV # 对目标进行统一最小化
    [levels, criLevel] = ea.ndsortDED(ObjV, None, 1, population.CV) # 只对个体划分出第一层
    [CombinObjV, weight] = ea.awGA(ObjV, population.CV) # 计算适应性权重以及多目标的加权单目标
    population.FitnV = (np.max(CombinObjV) - CombinObjV + 0.5) / (np.max(CombinObjV) - np.min(CombinObjV) + 0.5) # 计算种群适应度
    # 更新NDSet
    if NDSet is None:
        return population[np.where(levels == 1)[0]]
    else:
        tempPop = population[np.where(levels == 1)[0]] + NDSet # 将种群可行个体与NDSet合并
        [levels, criLevel] = ea.ndsortDED(maxormins * tempPop.ObjV, None, 1, tempPop.CV) # 只对个体划分出第一层
        liveIdx = np.where(levels == 1)[0] # 选择非支配个体
        NDSet = tempPop[liveIdx]
        # 对种群中被NDSet支配的个体进行惩罚
        punishFlag = np.zeros(population.sizes)
        punishFlag[np.where(liveIdx < population.sizes)[0]] == 1
        population.FitnV[np.where(punishFlag == 0)[0]] *= 0.5
        if len(liveIdx) > MAXSIZE: # 若要保留下来的NDSet个体数大于MAXSIZE，则根据拥挤距离进行筛选
            dis = ea.crowdis(NDSet.ObjV, levels[liveIdx]) # 计算拥挤距离
            NDSet = NDSet[ea.selecting('dup', np.array([dis]).T, MAXSIZE)] # 进行筛选
        return NDSet
    
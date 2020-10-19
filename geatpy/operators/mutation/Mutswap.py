# -*- coding: utf-8 -*-
from operators.mutation.Mutation import Mutation
from mutswap import mutswap


class Mutswap(Mutation):
    """
    Mutswap - class : 一个用于调用内核中的变异函数mutswap(染色体两点互换变异)的变异算子类，
                      该类的各成员属性与内核中的对应函数的同名参数含义一致，
                      可利用help(mutswap)查看各参数的详细含义及用法。
                      
    """

    def __init__(self, Pm=None, Parallel=False):
        self.Pm = Pm  # 表示染色体上变异算子所发生作用的最小片段发生变异的概率
        self.Parallel = Parallel  # 表示是否采用并行计算，缺省时默认为False

    def do(self, Encoding, OldChrom, FieldDR, *args):  # 执行变异
        return mutswap(Encoding, OldChrom, FieldDR, self.Pm, self.Parallel)

    def getHelp(self):  # 查看内核中的变异算子的API文档
        help(mutswap)

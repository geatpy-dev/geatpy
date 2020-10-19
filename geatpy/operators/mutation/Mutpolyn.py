# -*- coding: utf-8 -*-
from operators.mutation.Mutation import Mutation
from mutpolyn import mutpolyn


class Mutpolyn(Mutation):
    """
    Mutpolyn - class : 一个用于调用内核中的变异函数mutpolyn(多项式变异)的变异算子类，
                       该类的各成员属性与内核中的对应函数的同名参数含义一致，
                       可利用help(mutpolyn)查看各参数的详细含义及用法。
    """

    def __init__(self, Pm=None, DisI=20, FixType=1, Parallel=False):
        self.Pm = Pm  # 表示染色体上变异算子所发生作用的最小片段发生变异的概率
        self.DisI = DisI  # 多项式变异中的分布指数
        self.FixType = FixType  # 表示采用哪种方式来修复超出边界的染色体元素，可取值1，2，3，4，详细含义见help()帮助文档
        self.Parallel = Parallel  # 表示是否采用并行计算，缺省时默认为False

    def do(self, Encoding, OldChrom, FieldDR, *args):  # 执行变异
        return mutpolyn(Encoding, OldChrom, FieldDR, self.Pm, self.DisI, self.FixType, self.Parallel)

    def getHelp(self):  # 查看内核中的变异算子的API文档
        help(mutpolyn)

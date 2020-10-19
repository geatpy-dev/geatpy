# -*- coding: utf-8 -*-
from operators.mutation.Mutation import Mutation
from mutde import mutde


class Mutde(Mutation):
    """
    Mutde - class : 一个用于调用内核中的变异函数mutde(差分变异)的变异算子类，
                    该类的各成员属性与内核中的对应函数的同名参数含义一致，
                    可利用help(mutde)查看各参数的详细含义及用法。
                    
    """

    def __init__(self, F=0.5, FixType=1, Parallel=False):
        self.F = F  # 差分变异缩放因子
        self.FixType = FixType  # 表示采用哪种方式来修复超出边界的染色体元素，可取值1，2，3，4，详细含义见help()帮助文档
        self.Parallel = Parallel  # 表示是否采用并行计算，缺省时默认为False

    def do(self, Encoding, OldChrom, FieldDR, *args):  # 执行变异
        if len(args) != 0:
            XrList = args[0]
        else:
            XrList = None
        return mutde(Encoding, OldChrom, FieldDR, XrList, self.F, self.FixType, Parallel=self.Parallel)

    def getHelp(self):  # 查看内核中的变异算子的API文档
        help(mutde)

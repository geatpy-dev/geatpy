# -*- coding: utf-8 -*-
from operators.mutation.Mutation import Mutation
from mutpp import mutpp

class Mutpp(Mutation):
    """
    Mutpp - class : 一个用于调用内核中的变异函数mutpp(排列编码种群染色体变异)的变异算子类，
                    该类的各成员属性与内核中的对应函数的同名参数含义一致，
                    可利用help(mutpp)查看各参数的详细含义及用法。
    """
    def __init__(self, Pm = 1):
        self.Pm = Pm # 每条染色体发生变异的概率
    
    def do(self, Encoding, OldChrom, FieldDR, *args): # 执行变异
        return mutpp(Encoding, OldChrom, FieldDR, self.Pm)
    
    def getHelp(self): # 查看内核中的变异算子的API文档
        help(mutpp)
    
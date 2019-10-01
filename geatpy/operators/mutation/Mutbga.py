# -*- coding: utf-8 -*-
from operators.mutation.Mutation import Mutation
from mutbga import mutbga

class Mutbga(Mutation):
    """
    Mutbga - class : 一个用于调用内核中的变异函数mutbga(Breeder GA变异)的变异算子类，
                     该类的各成员属性与内核中的对应函数的同名参数含义一致，
                     可利用help(mutbga)查看各参数的详细含义及用法。
    """
    def __init__(self, Pm = 1, MutShrink = 0.5, Gradient = 20):
        self.Pm = Pm
        self.MutShrink = MutShrink
        self.Gradient = Gradient
    
    def do(self, Encoding, OldChrom, FieldDR, *args): # 执行变异
        return mutbga(Encoding, OldChrom, FieldDR, self.Pm, self.MutShrink, self.Gradient)
    
    def getHelp(self): # 查看内核中的变异算子的API文档
        help(mutbga)
    
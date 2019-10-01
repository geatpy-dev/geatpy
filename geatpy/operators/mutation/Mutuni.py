# -*- coding: utf-8 -*-
from operators.mutation.Mutation import Mutation
from mutuni import mutuni

class Mutuni(Mutation):
    """
    Mutuni - class : 一个用于调用内核中的变异函数mutuni(均匀变异)的变异算子类，
                     该类的各成员属性与内核中的对应函数的同名参数含义一致，
                     可利用help(mutuni)查看各参数的详细含义及用法。
    """
    def __init__(self, Pm = 1, Alpha = False, MutShrink = 1, Middle = False, Loop = False):
        self.Pm = Pm # 每条染色体发生变异的概率
        self.Alpha = Alpha # 表示均匀变异的变异半径
        self.MutShrink = MutShrink # 压缩率，用于压缩变异的范围
        self.Middle = Middle # 表示变异中心是否为搜索域的中央
        self.Loop = Loop # 表示是否采用循环的方式处理超出边界的变异结果
    
    def do(self, Encoding, OldChrom, FieldDR, *args): # 执行变异
        return mutuni(Encoding, OldChrom, FieldDR, self.Pm, self.Alpha, self.MutShrink, self.Middle, self.Loop)
    
    def getHelp(self): # 查看内核中的变异算子的API文档
        help(mutuni)
    
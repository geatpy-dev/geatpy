# -*- coding: utf-8 -*-
from operators.recombination.Recombination import Recombination
from recint import recint

class Recint(Recombination):
    
    """
    Recint - class : 一个用于调用内核中的函数recint(中间重组)的类，
                     该类的各成员属性与内核中的对应函数的同名参数含义一致，
                     可利用help(recint)查看各参数的详细含义及用法。
                     
    """
    
    def __init__(self, RecOpt = 0.7, Half = False, Parallel = False):
        self.RecOpt = RecOpt # 发生重组的概率
        self.Half = Half # 表示是否只保留一半重组结果
        self.Parallel = Parallel # 表示是否采用并行计算，缺省时默认为False
    
    def do(self, OldChrom): # 执行内核函数
        return recint(OldChrom, self.RecOpt, self.Half, self.Parallel)
    
    def getHelp(self): # 查看内核中的重组算子的API文档
        help(recint)
    
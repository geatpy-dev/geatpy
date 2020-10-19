# -*- coding: utf-8 -*-
from operators.recombination.Recombination import Recombination
from recndx import recndx


class Recndx(Recombination):
    """
    Recndx - class : 一个用于调用内核中的函数recndx(正态分布交叉)的类，
                     该类的各成员属性与内核中的对应函数的同名参数含义一致，
                     可利用help(recndx)查看各参数的详细含义及用法。
                     
    """

    def __init__(self, XOVR=0.7, Half_N=False, A=1.4826, Parallel=False):
        self.XOVR = XOVR  # 发生交叉的概率
        self.Half_N = Half_N  # 该参数是旧版的输入参数Half的升级版，用于控制交配方式
        self.A = A  # 一个跟交叉算子开发概率有关的参数
        self.Parallel = Parallel  # 表示是否采用并行计算，缺省时默认为False

    def do(self, OldChrom):  # 执行内核函数
        return recndx(OldChrom, self.XOVR, self.Half_N, self.A, self.Parallel)

    def getHelp(self):  # 查看内核中的重组算子的API文档
        help(recndx)

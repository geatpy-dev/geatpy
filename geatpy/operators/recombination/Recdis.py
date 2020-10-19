# -*- coding: utf-8 -*-
from operators.recombination.Recombination import Recombination
from recdis import recdis


class Recdis(Recombination):
    """
    Recdis - class : 一个用于调用内核中的函数recdis(离散重组)的类，
                     该类的各成员属性与内核中的对应函数的同名参数含义一致，
                     可利用help(recdis)查看各参数的详细含义及用法。
                     
    """

    def __init__(self, RecOpt=0.7, Half_N=False, GeneID=None, Parallel=False):
        self.RecOpt = RecOpt  # 发生重组的概率
        self.Half_N = Half_N  # 该参数是旧版的输入参数Half的升级版，用于控制交配方式
        self.GeneID = GeneID  # 基因ID，是一个行向量，若设置了该参数，则该函数会对具有相同基因ID的染色体片段进行整体离散重组
        self.Parallel = Parallel  # 表示是否采用并行计算，缺省时默认为False

    def do(self, OldChrom):  # 执行内核函数
        return recdis(OldChrom, self.RecOpt, self.Half_N, self.GeneID, self.Parallel)

    def getHelp(self):  # 查看内核中的重组算子的API文档
        help(recdis)

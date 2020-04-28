# -*- coding: utf-8 -*-
from operators.recombination.Recombination import Recombination
from xovmp import xovmp

class Xovmp(Recombination):
    
    """
    Xovmp - class : 一个用于调用内核中的函数xovmp(多点交叉)的类，
                    该类的各成员属性与内核中的对应函数的同名参数含义一致，
                    可利用help(xovmp)查看各参数的详细含义及用法。
                    
    """
    
    def __init__(self, XOVR = 0.7, Npt = 0, Half = False, GeneID = None):
        self.XOVR = XOVR # 发生交叉的概率
        self.Npt = Npt # 指明了交叉方式，Npt可为0(洗牌交叉)、1(单点交叉)、2(两点交叉)，若缺省或为None，则默认为0
        self.Half = Half # 表示是否只保留一半交叉结果
        self.GeneID = GeneID # 基因ID，是一个行向量，若设置了该参数，则该函数会对具有相同基因ID的染色体片段进行整体交叉。
    
    def do(self, OldChrom): # 执行内核函数
        return xovmp(OldChrom, self.XOVR, self.Npt, self.Half, self.GeneID)
    
    def getHelp(self): # 查看内核中的重组算子的API文档
        help(xovmp)
    
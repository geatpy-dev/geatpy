# -*- coding: utf-8 -*-

class Recombination:
    """
    Recombination - Interface : 重组算子接口
    进化算法框架中的所有重组算子类都实现该接口。
    
    """

    def __init__(self):
        pass

    def do(self):  # 用于调用内核层的重组算子函数执行重组
        pass

    def getHelp(self):  # 查看内核中的重组算子的API文档
        pass

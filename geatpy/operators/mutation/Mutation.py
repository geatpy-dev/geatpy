# -*- coding: utf-8 -*-

class Mutation:
    """
    Mutation - Interface : 变异算子接口
    进化算法框架中的所有变异算子类都实现该接口。
    
    """

    def __init__(self):
        pass

    def do(self):  # 用于调用内核中的变异算子函数执行变异
        pass

    def getHelp(self):  # 查看内核中的变异算子的API文档
        pass

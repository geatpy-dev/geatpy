# -*- coding: utf-8 -*-
from migrate import migrate


class Migrate:
    """
    Migrate - class : 一个用于调用内核中的种群迁移函数migrate的种群迁移算子类，
                      该类的各成员属性与内核中的对应函数的同名参数含义一致，
                      可利用help(migrate)查看各参数的详细含义及用法。
                      
    """

    def __init__(self, MIGR=0.2, Structure=0, Select=0, Replacement=0):
        self.MIGR = MIGR  # 表示种群个体的迁出概率
        self.Structure = Structure  # 表示种群迁移的结构（0：完全网状结构；1：邻近结构；2：环状结构）
        self.Select = Select  # 表示选择迁出个体的方式（0：随机选择迁出个体；1：择优选择迁出个体）
        self.Replacement = Replacement  # 表示种群在迁入个体时采用什么个体替换方式（0：替换迁出的个体；1：随机替换；2：择劣替换）

    def do(self, populations, *args):  # 执行变异，populations为存储着种群类对象的列表
        if type(populations) != list:
            raise RuntimeError('error in Migrate: The populations must be a list. (输入参数populations必须是list类型。)')

        PopSizes = list(pop.sizes for pop in populations)
        FitnVs = list(pop.FitnV for pop in populations)
        # 调用种群迁移算子进行种群个体迁移
        [Aborigines, Foreigners, FromPlaces] = migrate(PopSizes, self.MIGR, self.Structure, self.Select,
                                                       self.Replacement, FitnVs)
        NewPopulations = []
        for i in range(len(populations)):  # 更新迁移个体后的种群
            NewPopulations.append((populations[i])[Aborigines[i]] + (populations[FromPlaces[i]])[Foreigners[i]])
        return NewPopulations

    def getHelp(self):  # 查看内核中的变异算子的API文档
        help(migrate)

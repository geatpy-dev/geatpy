# -*- coding: utf-8 -*-
import os
import numpy as np
import geatpy as ea

class PsyPopulation(ea.Population):
    
    """
PsyPopulation : class - 多染色体种群类(Popysomy Population)

描述:
    多染色体种群类是用来存储每个个体包含多条染色体的种群相关信息的一个类。
    该类和种群类Population似，不同之处是可以包含多条染色体，因此支持复杂的混合编码。
    
属性:
    sizes     : int   - 种群规模，即种群的个体数目。
    
    ChromNum  : int   - 染色体的数目，即每个个体有多少条染色体。
    
    Linds     : list  - 存储种群各染色体长度的列表。
    
    Encodings : list  - 存储各染色体编码方式的列表。
    
    Fields    : list  - 存储各染色体对应的译码矩阵的列表。
    
    Chroms    : list  - 存储种群各染色体矩阵的列表。
    
    ObjV      : array - 种群目标函数值矩阵，每一行对应一个个体的目标函数值，每一列对应一个目标。
    
    FitnV     : array - 种群个体适应度列向量，每个元素对应一个个体的适应度，最小适应度为0。
    
    CV        : array - CV(Constraint Violation Value)是用来定量描述违反约束条件程度的矩阵，每行对应一个个体，每列对应一个约束。
                        注意：当没有用到约束条件时，种群也会携带一个只有一列的、元素全为0的CV。
    
    Phen      : array - 种群表现型矩阵（即染色体解码后所代表的决策变量所组成的矩阵）。
    
函数:
    详见源码。

"""

    def __init__(self, Encodings, Fields, NIND, Chroms = None, ObjV = None, FitnV = None, CV = None, Phen = None):
        """
        描述: 种群类的构造方法，用于实例化种群对象，例如：
             import geatpy as ea
             population = ea.PsyPopulation(Encodings, Fields, NIND)，
             NIND为所需要的个体数，
             此时得到的population还没被真正初始化，仅仅是完成种群对象的实例化。
             该构造方法必须传入Chroms，才算是完成种群真正的初始化。
             一开始可以只传入Encodings, Fields以及NIND来完成种群对象的实例化，
             其他属性可以后面再通过计算进行赋值。
        """
        
        self.sizes = NIND
        self.ChromNum = len(Encodings)
        self.Linds = []
        self.Encodings = Encodings
        self.Fields = Fields.copy()
        self.Chroms = [None] * self.ChromNum # 初始化Chroms为元素是None的列表
        if Chroms is None:
            self.Linds = [0] * self.ChromNum
        else:
            for i in range(self.ChromNum):
                if Chroms[i] is not None:
                    self.Linds.append(Chroms[i].shape[1])
                    self.Chroms[i] = Chroms[i].copy() if Chroms[i] is not None else Chroms[i]
                else:
                    self.Linds.append(0)
        self.ObjV = ObjV.copy() if ObjV is not None else ObjV
        self.FitnV = FitnV.copy() if FitnV is not None else np.ones((self.sizes, 1))
        self.CV = CV.copy() if CV is not None else np.zeros((self.sizes, 1))
        self.Phen = Phen.copy() if Phen is not None else Phen
    
    def initChrom(self, NIND = None):
        """
        描述: 初始化种群染色体矩阵，NIND为所需要的个体数，并返回一个新的种群对象。
        NIND可缺省，不缺省时，种群在初始化染色体矩阵前会把种群规模调整为NIND。
        """
        
        if NIND is not None:
            self.sizes = NIND # 重新设置种群规模
        # 遍历各染色体矩阵进行初始化
        for i in range(self.ChromNum):
            self.Chroms[i] = ea.crtpc(self.Encodings[i], self.sizes, self.Fields[i]) # 生成染色体矩阵
        self.ObjV = None
        self.FitnV = np.ones((self.sizes, 1))
        self.CV = np.zeros((self.sizes, 1))
        self.Phen = self.decoding() # 解码
    
    def decoding(self):
        """
        描述: 种群染色体解码。
        """
        
        Phen = np.ones((self.sizes, 0)) # 初始化一个空的矩阵
        # 遍历各染色体矩阵进行解码
        for i in range(self.ChromNum):
            if self.Encodings[i] == 'BG': # 此时Field实际上为FieldD
                tempPhen = ea.bs2ri(self.Chroms[i], self.Fields[i]) # 把二进制转化为实值
            elif self.Encodings[i] == 'RI' or self.Encodings[i] == 'P':
                tempPhen = self.Chroms[i].copy()
            else:
                raise RuntimeError('error in Population.decoding: Encoding must be ''BG'' or ''RI'' or ''P''. (编码设置有误，Encoding必须为''BG'', ''RI'' 或 ''P''。)')
            Phen = np.hstack([Phen, tempPhen])
        return Phen
    
    def copy(self):
        """
        copy : function - 种群的复制
        用法:
            假设pop是一个种群矩阵，那么：pop1 = pop.copy()即可完成对pop种群的复制。
        """
        
        return PsyPopulation(self.Encodings, 
                          self.Fields, 
                          self.sizes, 
                          self.Chroms, 
                          self.ObjV, 
                          self.FitnV, 
                          self.CV, 
                          self.Phen)
    
    def __getitem__(self, index):
        """
        描述: 种群的切片，即根据index下标向量选出种群中相应的个体组成一个新的种群。
        用法: 假设pop是一个包含多于2个个体的种群矩阵，那么：
             pop1 = pop[[0,1]]即可得到由pop种群的第1、2个个体组成的种群。
        注意: index必须是一个Numpy array类型的行向量。
        """
        
        NewChroms = []
        for i in range(self.ChromNum):
            if self.Chroms[i] is None:
                raise RuntimeError('error in Population: Chrom is None. (种群染色体矩阵未初始化。)')
            NewChroms.append(self.Chroms[i][index])
        NIND = NewChroms[0].shape[0]
        return PsyPopulation(self.Encodings, 
                             self.Fields, 
                             NIND,
                             NewChroms, 
                             self.ObjV[index], 
                             self.FitnV[index], 
                             self.CV[index], 
                             self.Phen[index])
    
    def shuffle(self):
        """
        shuffle : function - 打乱种群个体的个体顺序
        用法: 假设pop是一个种群矩阵，那么，pop.shuffle()即可完成对pop种群个体顺序的打乱
        """
        
        shuff = np.argsort(np.random.rand(self.sizes))
        for i in range(self.ChromNum):
            if self.Chroms[i] is None:
                raise RuntimeError('error in Population: Chrom is None. (种群染色体矩阵未初始化。)')
            self.Chroms[i] = self.Chroms[i][shuff, :]
        self.ObjV = self.ObjV[shuff, :] if self.ObjV is not None else self.ObjV
        self.FitnV = self.FitnV[shuff]
        self.CV = self.CV[shuff, :]
        self.Phen = self.Phen[shuff, :]
    
    def __setitem__(self, index, pop): # 种群个体赋值
        """
        描述: 种群个体的赋值
        用法: 假设pop是一个包含多于2个个体的种群矩阵，pop1是另一个包含2个个体的种群矩阵，那么
             pop[[0,1]] = pop1，即可完成将pop种群的第1、2个个体赋值为pop1种群的个体。
        注意: index必须是一个Numpy array类型的行向量。
        """
        
        if self.sizes != pop.sizes:
            raise RuntimeError('error in Population: Sizes disagree. (两者的规模必须一致。)')
        for i in range(self.ChromNum):
            if self.Encodings[i] != pop.Encodings[i]:
                raise RuntimeError('error in Population: Encoding disagree. (两种群染色体的编码方式必须一致。)')
            if np.all(self.Fields[i] == pop.Fields[i]) == False:
                raise RuntimeError('error in Population: Field disagree. (两者的译码矩阵必须一致。)')
            if self.Chroms[i] is None:
                raise RuntimeError('error in Population: Chrom is None. (种群染色体矩阵未初始化。)')
            self.Chroms[i][index] = pop.Chroms[i]
        self.sizes = self.Chrom.shape[0] # 更新种群规模
        self.ObjV[index] = pop.ObjV
        self.FitnV = np.ones((self.sizes, 1)) # 重置适应度
        self.CV[index] = pop.CV
        self.Phen[index] = pop.Phen
    
    def __add__(self, pop):
        """
        描述: 种群个体合并
        用法: 假设pop1, pop2是两个种群，它们的个体数可以相等也可以不相等，此时
             pop = pop1 + pop2，即可完成对pop1和pop2两个种群个体的合并
        """
        
        NIND = self.sizes + pop.sizes # 得到合并种群的个体数
        NewChroms = self.Chroms
        for i in range(self.ChromNum):
            if self.Encodings[i] != pop.Encodings[i]:
                raise RuntimeError('error in Population: Encoding disagree. (两种群染色体的编码方式必须一致。)')
            if np.all(self.Fields[i] == pop.Fields[i]) == False:
                raise RuntimeError('error in Population: Field disagree. (两者的译码矩阵必须一致。)')
            if self.Chroms[i] is None or pop.Chroms[i] is None:
                raise RuntimeError('error in Population: Chrom is None. (种群染色体矩阵未初始化。)')
            NewChroms[i] = np.vstack([NewChroms[i], pop.Chroms[i]])
        return PsyPopulation(self.Encodings, 
                             self.Fields, 
                             NIND, 
                             NewChroms, 
                             np.vstack([self.ObjV, pop.ObjV]), 
                             np.ones((NIND, 1)), # 重置适应度
                             np.vstack([self.CV, pop.CV]), 
                             np.vstack([self.Phen, pop.Phen]))
    
    def __len__(self):
        """
        描述: 计算种群规模
        用法: 假设pop是一个种群，那么len(pop)即可得到该种群的个体数。
             实际上，种群规模也可以通过pop.sizes得到
        """
        
        return self.sizes
    
    def save(self):
        """
        描述: 把种群的信息保存到文件中。
        该函数将在"Result"文件夹下保存种群的信息，其中：
        "Encodingsi.txt"保存种群的染色体编码，i为0,1,2,3...；
        "Fieldsi.csv"保存种群染色体的译码矩阵，i为0,1,2,3...；
        "Chromsi.csv"保存种群的染色体矩阵，i为0,1,2,3...；
        "ObjV.csv"保存种群的目标函数矩阵；
        "FitnV.csv"保存种群个体的适应度列向量；
        "CV.csv"保存种群个体的违反约束程度矩阵；
        "Phen.csv"保存种群染色体表现型矩阵；
        注意：该函数不会对种群的合法性进行检查。
        """
        
        if os.path.exists('Result') == False:
            os.makedirs('Result')
        for i in range(self.ChromNum):
            with open('Result/Encodings' + str(i) + '.txt','w') as file:
                file.write(self.Encodings[i])
            np.savetxt('Result/Fields' + str(i) + '.csv', self.Fields[i], delimiter=',')
            np.savetxt('Result/Chroms' + str(i) + '.csv', self.Chroms[i], delimiter=',')
        np.savetxt('Result/ObjV.csv', self.ObjV, delimiter=',')
        np.savetxt('Result/FitnV.csv', self.FitnV, delimiter=',')
        np.savetxt('Result/CV.csv', self.CV, delimiter=',')
        np.savetxt('Result/Phen.csv', self.Phen, delimiter=',')
        print('种群信息导出完毕。')
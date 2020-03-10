# -*- coding: utf-8 -*-
import os
import numpy as np
import geatpy as ea

class Population:
    
    """
Population : class - 种群类

描述:
    种群类是用来存储种群相关信息的一个类。

属性:
    sizes    : int   - 种群规模，即种群的个体数目。
    
    Encoding : str   - 染色体编码方式，
                       'BG':二进制/格雷编码；
                       'RI':实整数编码，即实数和整数的混合编码；
                       'P':排列编码。
                       相关概念：术语“实值编码”包含实整数编码和排列编码，
                       它们共同的特点是染色体不需要解码即可直接表示对应的决策变量。
                       "实整数"指的是种群染色体既包含实数的小数，也包含实数的整数。
                       特殊用法：
                       设置Encoding=None，此时种群类的Field,Chrom成员属性将被设置为None，
                       种群将不携带与染色体直接相关的信息，可以减少不必要的数据存储，
                       这种用法可以在只想统计非染色体直接相关的信息时使用，
                       尤其可以在多种群进化优化过程中对个体进行统一的适应度评价时使用。
    
    Field    : array - 译码矩阵，可以是FieldD或FieldDR（详见Geatpy数据结构）。
    
    Chrom    : array - 种群染色体矩阵，每一行对应一个个体的一条染色体。
    
    Lind     : int   - 种群染色体长度。
    
    ObjV     : array - 种群目标函数值矩阵，每一行对应一个个体的目标函数值，每一列对应一个目标。
    
    FitnV    : array - 种群个体适应度列向量，每个元素对应一个个体的适应度，最小适应度为0。
    
    CV       : array - CV(Constraint Violation Value)是用来定量描述违反约束条件程度的矩阵，每行对应一个个体，每列对应一个约束。
                       注意：当没有用到约束条件时，种群也会携带一个只有一列的、元素全为0的CV。
    
    Phen     : array - 种群表现型矩阵（即种群各染色体解码后所代表的决策变量所组成的矩阵）。
    
函数:
    详见源码。

"""

    def __init__(self, Encoding, Field, NIND, Chrom = None, ObjV = None, FitnV = None, CV = None, Phen = None):
        """
        描述: 种群类的构造方法，用于实例化种群对象，例如：
             import geatpy as ea
             population = ea.Population(Encoding, Field, NIND)，
             NIND为所需要的个体数，
             此时得到的population还没被真正初始化，仅仅是完成种群对象的实例化。
             该构造方法必须传入Chrom，才算是完成种群真正的初始化。
             一开始可以只传入Encoding, Field以及NIND来完成种群对象的实例化，
             其他属性可以后面再通过计算进行赋值。
        """
        
        if type(NIND) is int and NIND >= 0:
            self.sizes = NIND 
        else:
            raise RuntimeError('error in Population: Size error. (种群规模设置有误，必须为非负整数。)')
        self.Encoding = Encoding
        if Encoding is None:
            self.Field = None
            self.Chrom = None
        else:
            self.Field = Field.copy()
            self.Chrom = Chrom.copy() if Chrom is not None else Chrom
        self.Lind = Chrom.shape[1] if Chrom is not None else 0
        self.ObjV = ObjV.copy() if ObjV is not None else ObjV
        self.FitnV = FitnV.copy() if FitnV is not None else np.ones((self.sizes, 1))
        self.CV = CV.copy() if CV is not None else np.zeros((self.sizes, 1))
        self.Phen = Phen.copy() if Phen is not None else Phen
    
    def initChrom(self, NIND = None):
        """
        描述: 初始化种群染色体矩阵，NIND为所需要的个体数。
        NIND可缺省，不缺省时，种群在初始化染色体矩阵前会把种群规模调整为NIND。
        """
        
        if NIND is not None:
            self.sizes = NIND # 重新设置种群规模
        self.Chrom = ea.crtpc(self.Encoding, self.sizes, self.Field) # 生成染色体矩阵
        self.Lind = self.Chrom.shape[1] # 计算染色体的长度
        self.ObjV = None
        self.FitnV = np.ones((self.sizes, 1))
        self.CV = np.zeros((self.sizes, 1))
        self.Phen = self.decoding() # 解码
    
    def setChrom(self, Chrom = None):
        """
        描述: 该函数用于插入种群染色体矩阵。
        当想要插入一些染色体时，可以调用该函数（例如想要插入一些已知的较优解对应的染色体）。
        调用该函数之前，种群染色体必须已经被初始化。
        该函数会简单地把现有的种群染色体全部或部分替换成待插入的染色体，不会择劣替换。
        """
        
        print('Warning: Population.setChrom() will be removed in the next version. (警告：Population.setChrom()函数将会在后面的版本中被移除。若想要在进化前插入先验知识，请使用Algorithm类新增的“先知种群”。)')
        
        if Chrom is None: # 如果Chrom是None，表示不需要插入，直接退出
            return
        if self.Chrom is None:
            raise RuntimeError('error in Population.setChrom: Chrom is None. (种群染色体矩阵未初始化。)')
        NewChrom = np.vstack([Chrom, self.Chrom])
        self.Chrom = NewChrom[:self.sizes, :]
        
    def decoding(self):
        """
        描述: 种群染色体解码。
        """
    
        if self.Encoding == 'BG': # 此时Field实际上为FieldD
            Phen = ea.bs2ri(self.Chrom, self.Field) # 把二进制转化为实整数
        elif self.Encoding == 'RI' or self.Encoding == 'P':
            Phen = self.Chrom.copy()
        else:
            raise RuntimeError('error in Population.decoding: Encoding must be ''BG'' or ''RI'' or ''P''. (编码设置有误，解码时Encoding必须为''BG'', ''RI'' 或 ''P''。)')
        return Phen
    
    def copy(self):
        """
        copy : function - 种群的复制
        用法:
            假设pop是一个种群矩阵，那么：pop1 = pop.copy()即可完成对pop种群的复制。
        """
        
        return Population(self.Encoding, 
                          self.Field, 
                          self.sizes, 
                          self.Chrom, 
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
        
        if self.Encoding is None:
            NewChrom = None
        else:
            if self.Chrom is None:
                raise RuntimeError('error in Population: Chrom is None. (种群染色体矩阵未初始化。)')
            NewChrom = self.Chrom[index]
        NewFitnV = self.FitnV[index]
        NIND = NewFitnV.shape[0]
        return Population(self.Encoding, 
                          self.Field, 
                          NIND,
                          NewChrom, 
                          self.ObjV[index] if self.ObjV is not None else self.ObjV, 
                          NewFitnV, 
                          self.CV[index], 
                          self.Phen[index])
    
    def shuffle(self):
        """
        shuffle : function - 打乱种群个体的个体顺序
        用法: 假设pop是一个种群矩阵，那么，pop.shuffle()即可完成对pop种群个体顺序的打乱
        """
        
        shuff = np.argsort(np.random.rand(self.sizes))
        if self.Encoding is None:
            self.Chrom = None
        else:
            if self.Chrom is None:
                raise RuntimeError('error in Population: Chrom is None. (种群染色体矩阵未初始化。)')
            self.Chrom = self.Chrom[shuff, :]
        self.ObjV = self.ObjV[shuff, :] if self.ObjV is not None else self.ObjV
        self.FitnV = self.FitnV[shuff]
        self.CV = self.CV[shuff, :]
        self.Phen = self.Phen[shuff, :]
    
    def __setitem__(self, index, pop): # 种群个体赋值（种群个体替换）
        """
        描述: 种群个体的赋值
        用法: 假设pop是一个包含多于2个个体的种群矩阵，pop1是另一个包含2个个体的种群矩阵，那么
             pop[[0,1]] = pop1，即可完成将pop种群的第1、2个个体赋值为pop1种群的个体。
        注意: index必须是一个Numpy array类型的行向量。
             此外，进行种群个体替换后，该函数不会对适应度进行主动重置，
             如果因个体替换而需要重新对所有个体的适应度进行评价，则需要手写代码更新种群的适应度。
        """
        
        if self.Encoding is not None:
            if self.Encoding != pop.Encoding:
                raise RuntimeError('error in Population: Encoding disagree. (两种群染色体的编码方式必须一致。)')
            if np.all(self.Field == pop.Field) == False:
                raise RuntimeError('error in Population: Field disagree. (两者的译码矩阵必须一致。)')
            if len(index) != pop.sizes:
                raise RuntimeError('error in Population: Sizes disagree. (index的长度和待插入的种群规模必须一致。)')
            if self.Chrom is None:
                raise RuntimeError('error in Population: Chrom is None. (种群染色体矩阵未初始化。)')
            self.Chrom[index] = pop.Chrom
        self.ObjV[index] = pop.ObjV
        self.FitnV[index] = pop.FitnV
        self.CV[index] = pop.CV
        self.Phen[index] = pop.Phen
        self.sizes = self.Phen.shape[0] # 更新种群规模
    
    def __add__(self, pop):
        """
        描述: 种群个体合并
        用法: 假设pop1, pop2是两个种群，它们的个体数可以相等也可以不相等，此时
             pop = pop1 + pop2，即可完成对pop1和pop2两个种群个体的合并。
        注意：
            进行种群合并后，该函数不会对适应度进行主动重置，
            如果因种群合并而需要重新对所有个体的适应度进行评价，则需要手写代码更新种群的适应度。
        """
        
        if self.Encoding is None:
            NewChrom = None
        else:
            if self.Encoding != pop.Encoding:
                raise RuntimeError('error in Population: Encoding disagree. (两种群染色体的编码方式必须一致。)')
            if np.all(self.Field == pop.Field) == False:
                raise RuntimeError('error in Population: Field disagree. (两者的译码矩阵必须一致。)')
            if self.Chrom is None or pop.Chrom is None:
                raise RuntimeError('error in Population: Chrom is None. (种群染色体矩阵未初始化。)')
            NewChrom = np.vstack([self.Chrom, pop.Chrom])
        NIND = self.sizes + pop.sizes # 得到合并种群的个体数
        return Population(self.Encoding, 
                          self.Field, 
                          NIND, 
                          NewChrom, 
                          np.vstack([self.ObjV, pop.ObjV]), 
                          np.vstack([self.FitnV, pop.FitnV]),
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
        "Encoding.txt"保存种群的染色体编码；
        "Field.csv"保存种群染色体的译码矩阵；
        "Chrom.csv"保存种群的染色体矩阵；
        "ObjV.csv"保存种群的目标函数矩阵；
        "FitnV.csv"保存种群个体的适应度列向量；
        "CV.csv"保存种群个体的违反约束程度矩阵；
        "Phen.csv"保存种群染色体表现型矩阵；
        注意：该函数不会对种群的合法性进行检查。
        """
        
        if os.path.exists('Result') == False:
            os.makedirs('Result')
        with open('Result/Encoding.txt','w') as file:
            file.write(str(self.Encoding))
        if self.Encoding is not None:
            np.savetxt('Result/Field.csv', self.Field, delimiter=',')
            np.savetxt('Result/Chrom.csv', self.Chrom, delimiter=',')
        np.savetxt('Result/ObjV.csv', self.ObjV, delimiter=',')
        np.savetxt('Result/FitnV.csv', self.FitnV, delimiter=',')
        np.savetxt('Result/CV.csv', self.CV, delimiter=',')
        np.savetxt('Result/Phen.csv', self.Phen, delimiter=',')
        print('种群信息导出完毕。')
    
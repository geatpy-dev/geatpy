# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ga

class Population:
    
    """
Population : class - 种群类

描述:
    种群类是用来存储种群相关信息的一个类，注意：种群的所有个体的染色体的编码方式都是相同的，
    若需要采用混合编码，需要使用多个种群对象。

属性:
    sizes    : int   - 种群规模，即种群的个体数目。
    
    Lind     : int   - 种群染色体长度。
    
    Encoding : str   - 染色体编码方式，
                       'B':二进制编码; 'G':格雷编码; 'R':实数编码; 'I':整数编码; 'P':排列编码。
                     
    conordis : int   - 连续或离散标记，0表示该种群染色体解码后得到的决策变量是连续的；
                       1表示该种群染色体解码后得到的变量是离散的。
    
    Field    : array - 译码矩阵，可以是FieldD或FieldDR（详见Geatpy数据结构）。
    
    Chrom    : array - 种群染色体矩阵，每一行对应一个个体的一条染色体。
    
    ObjV     : array - 种群目标函数值矩阵，每一行对应一个个体的目标函数值，每一列对应一个目标。
    
    FitnV    : array - 种群个体适应度列向量，每个元素对应一个个体的适应度，最小适应度为0。
    
    CV       : array - CV(Constraint Violation Value)是用来定量描述违反约束条件程度的矩阵，每行对应一个个体，每列对应一个约束。
    
    Phen     : array - 种群表现型矩阵（即种群各染色体解码后所代表的决策变量所组成的矩阵）。
    
函数:
    详见源码

"""

    def __init__(self, Encoding, conordis, Field, NIND, Chrom = None, ObjV = None, FitnV = None, CV = None, Phen = None):
        """
        描述: 种群类的构造方法，用于实例化种群对象，例如：population = ga.Population(Encoding, conordis, Field, NIND)，
             此时得到的population还没被真正初始化，只是设置了种群的一些基本属性，
             该构造方法必须传入Chrom，才算是完成种群真正的初始化。
             NIND为所需要的个体数，
             一开始可以只传入Encoding, conordis, Field以及NUM，此时只设置了种群的一些基本属性，
             其他属性可以后面再通过其他计算进行赋值
        """
        
        self.sizes = NIND
        self.Lind = Chrom.shape[1] if Chrom is not None else 0
        self.Encoding = Encoding
        self.conordis = conordis
        self.Field = Field.copy()
        self.Chrom = Chrom.copy() if Chrom is not None else Chrom
        self.ObjV = ObjV.copy() if ObjV is not None else ObjV
        if FitnV is None: # 若没传入适应度，则个体的适应度默认为1
            self.FitnV = np.ones((self.sizes, 1))
        else:
            self.FitnV = FitnV.copy()
        if CV is None: # 若没传入CV，则个体的CV默认为0矩阵
            CV = np.zeros((self.sizes, 1))
        else:
            self.CV = CV.copy()
        self.Phen = Phen.copy() if Phen is not None else Phen
    
    def initChrom(self, NIND):
        """
        描述: 初始化种群染色体矩阵，NIND为所需要的个体数，并返回一个新的种群对象
        """
        
        Lind = 0
        self.sizes = NIND
        if self.Encoding == 'B' or self.Encoding == 'G': # 生成二进制/格雷码种群，此时Field实际上为FieldD
            Lind = int(np.sum(self.Field[0, :])) # 根据Field种群染色体长度
            self.Chrom = ga.crtip(self.sizes, np.vstack([np.zeros((1, Lind)), np.ones((1, Lind))]))
        elif self.Encoding == 'R':
            self.Chrom = ga.crtrp(self.sizes, self.Field) # 生成实数值种群，此时Field实际上为FieldDR
        elif self.Encoding == 'I':
            self.Chrom = ga.crtip(self.sizes, self.Field) # 生成整数值种群，此时Field实际上为FieldDR
        elif self.Encoding == 'P':
            self.Chrom = ga.crtpp(self.sizes, self.Field) # 生成排列编码种群
        else:
            raise RuntimeError('error in Population.initChrom: Encoding is illegal. (Encoding参数错误。)')
        self.FitnV = np.ones((self.sizes, 1))
        self.CV = np.zeros((self.sizes, 1))
        self.Phen = self.decoding() # 解码
    
    def decoding(self):
        """
        描述: 种群染色体解码。
    """
    
        if self.Encoding == 'B' or self.Encoding == 'G': # 此时Field实际上为FieldD
            if self.conordis == 0:
                Phen = ga.bs2rv(self.Chrom, self.Field) # 解码
            elif self.conordis == 1:
                Phen = ga.bs2int(self.Chrom, self.Field)
            else:
                raise RuntimeError('error in Population.decoding: Conordis is illegal. (种群的conordis参数错误。)')
        else:
            Phen = self.Chrom.copy()
        return Phen
    
    def copy(self):
        """
        copy : function - 种群的复制
        用法:
            假设pop是一个种群矩阵，那么：pop1 = pop.copy()即可完成对pop种群的复制。
        """
        
        return Population(self.Encoding, self.conordis, self.Field, self.sizes, self.Chrom, self.ObjV, self.FitnV, self.CV, self.Phen)
    
    def __getitem__(self, index):
        """
        描述: 种群的切片，即根据index下标向量选出种群中相应的个体组成一个新的种群。
        用法: 假设pop是一个包含多于2个个体的种群矩阵，那么：
             pop1 = pop[[0,1]]即可得到由pop种群的第1、2个个体组成的种群。
        注意: index必须是一个Numpy array类型的行向量。
        """
        
        NIND = len(index)
        if self.sizes < NIND:
            raise RuntimeError('error in Population: Index is out of bounds. (索引越界。)')
        if self.Chrom is None:
            raise RuntimeError('error in Population: Chrom is None. (种群染色体矩阵未初始化。)')
        return Population(self.Encoding, self.conordis, self.Field, NIND, self.Chrom[index], self.ObjV[index], self.FitnV[index], self.CV[index], self.Phen[index])
    
    def shuffle(self):
        """
        shuffle : function - 打乱种群个体的个体顺序
        用法: 假设pop是一个种群矩阵，那么，pop.shuffle()即可完成对pop种群个体顺序的打乱
        """
        
        shuff = np.argsort(np.random.rand(self.sizes))
        self.Chrom = self.Chrom[shuff, :]
        if self.ObjV is not None:
            self.ObjV = self.ObjV[shuff, :]
        self.FitnV = self.FitnV[shuff]
        self.CV = self.CV[shuff, :]
        self.Phen = self.Phen[shuff, :]
    
    def __setitem__(self, index, pop): # 种群个体赋值
        """
        描述: 种群个体的赋值
        用法: 假设pop是一个包含多于2个个体的种群矩阵，pop1是另一个包含2个个体的种群矩阵，那么
             pop[[0,1]] = pop1，即可完成将pop种群的第1、2个个体赋值为pop1种群的个体。
        """
        
        if self.Encoding != pop.Encoding:
            raise RuntimeError('error in Population: Encoding disagree. (两种群染色体的编码方式必须一致。)')
        if self.conordis != pop.conordis:
            raise RuntimeError('error in Population: Conordis disagree. (两种群染色体所代表的变量的连续或离散性必须一致。)')
        if np.all(self.Field == pop.Field) == False:
            raise RuntimeError('error in Population: Field disagree. (两者的译码矩阵必须一致。)')
        if self.sizes != pop.sizes:
            raise RuntimeError('error in Population: Sizes disagree. (两者的规模必须一致。)')
        self.Chrom[index] = pop.Chrom
        self.ObjV[index] = pop.ObjV
        self.FitnV[index] = pop.FitnV
        self.CV[index] = pop.CV
        self.Phen[index] = pop.Phen
    
    def __add__(self, pop):
        """
        描述: 种群个体合并
        用法: 假设pop1, pop2是两个种群，它们的个体数可以相等也可以不相等，此时
             pop = pop1 + pop2，即可完成对pop1和pop2两个种群个体的合并
        """
        
        if self.Encoding != pop.Encoding:
            raise RuntimeError('error in Population: Encoding disagree. (两种群染色体的编码方式必须一致。)')
        if self.conordis != pop.conordis:
            raise RuntimeError('error in Population: Conordis disagree. (两种群染色体所代表的变量的连续或离散性必须一致。)')
        if np.all(self.Field == pop.Field) == False:
            raise RuntimeError('error in Population: Field disagree. (两者的译码矩阵必须一致。)')
        return Population(self.Encoding, self.conordis, self.Field, self.sizes + pop.sizes, np.vstack([self.Chrom, pop.Chrom]), np.vstack([self.ObjV, pop.ObjV]), np.vstack([self.FitnV, pop.FitnV]), np.vstack([self.CV, pop.CV]), np.vstack([self.Phen, pop.Phen]))

    def __len__(self):
        """
        描述: 计算种群规模
        用法: 假设pop是一个种群，那么len(pop)即可得到该种群的个体数。
             实际上，种群规模也可以通过pop.sizes得到
        """
        
        return self.sizes
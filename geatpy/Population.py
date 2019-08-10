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
    
    Lind     : int   - 种群染色体长度。
    
    Encoding : str   - 染色体编码方式，
                       'BG':二进制/格雷编码；
                       'RI':实整数编码，即实数和整数的混合编码；
                       'P':排列编码。
                       相关概念：术语“实值编码”包含实整数编码和排列编码，
                       它们共同的特点是染色体不需要解码即可直接表示对应的决策变量。
                       "实整数"指的是种群染色体既包含实数的小数，也包含实数的整数。
    
    Field    : array - 译码矩阵，可以是FieldD或FieldDR（详见Geatpy数据结构）。
    
    Chrom    : array - 种群染色体矩阵，每一行对应一个个体的一条染色体。
    
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
        
        self.sizes = NIND
        self.Lind = Chrom.shape[1] if Chrom is not None else 0
        self.Encoding = Encoding
        self.Field = Field.copy()
        self.Chrom = Chrom.copy() if Chrom is not None else Chrom
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
        self.Chrom = ea.crtpc(self.Encoding, self.sizes, self.Field) # 生成染色体矩阵
        self.ObjV = None
        self.FitnV = np.ones((self.sizes, 1))
        self.CV = np.zeros((self.sizes, 1))
        self.Phen = self.decoding() # 解码
    
    def decoding(self):
        """
        描述: 种群染色体解码。
        """
    
        if self.Encoding == 'BG': # 此时Field实际上为FieldD
            Phen = ea.bs2ri(self.Chrom, self.Field) # 把二进制转化为实值
        elif self.Encoding == 'RI' or self.Encoding == 'P':
            Phen = self.Chrom.copy()
        else:
            raise RuntimeError('error in Population.decoding: Encoding must be ''BG'' or ''RI'' or ''P''. (编码设置有误，Encoding必须为''BG'', ''RI'' 或 ''P''。)')
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
        
        if self.Chrom is None:
            raise RuntimeError('error in Population: Chrom is None. (种群染色体矩阵未初始化。)')
        NewChrom = self.Chrom[index]
        NIND = NewChrom.shape[0]
        return Population(self.Encoding, 
                          self.Field, 
                          NIND,
                          NewChrom, 
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
        
        if self.Chrom is None:
            raise RuntimeError('error in Population: Chrom is None. (种群染色体矩阵未初始化。)')
        self.Chrom = self.Chrom[shuff, :]
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
        
        if self.Encoding != pop.Encoding:
            raise RuntimeError('error in Population: Encoding disagree. (两种群染色体的编码方式必须一致。)')
        if np.all(self.Field == pop.Field) == False:
            raise RuntimeError('error in Population: Field disagree. (两者的译码矩阵必须一致。)')
        if self.sizes != pop.sizes:
            raise RuntimeError('error in Population: Sizes disagree. (两者的规模必须一致。)')
        if self.Chrom is None:
            raise RuntimeError('error in Population: Chrom is None. (种群染色体矩阵未初始化。)')
        self.Chrom[index] = pop.Chrom
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
        
        if self.Encoding != pop.Encoding:
            raise RuntimeError('error in Population: Encoding disagree. (两种群染色体的编码方式必须一致。)')
        if np.all(self.Field == pop.Field) == False:
            raise RuntimeError('error in Population: Field disagree. (两者的译码矩阵必须一致。)')
        if self.Chrom is None or pop.Chrom is None:
            raise RuntimeError('error in Population: Chrom is None. (种群染色体矩阵未初始化。)')
        NIND = self.sizes + pop.sizes # 得到合并种群的个体数
        return Population(self.Encoding, 
                          self.Field, 
                          NIND, 
                          np.vstack([self.Chrom, pop.Chrom]), 
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
            file.write(self.Encoding)
        np.savetxt('Result/Field.csv', self.Field, delimiter=',')
        np.savetxt('Result/Chrom.csv', self.Chrom, delimiter=',')
        np.savetxt('Result/ObjV.csv', self.ObjV, delimiter=',')
        np.savetxt('Result/FitnV.csv', self.FitnV, delimiter=',')
        np.savetxt('Result/CV.csv', self.CV, delimiter=',')
        np.savetxt('Result/Phen.csv', self.Phen, delimiter=',')
        print('种群信息导出完毕。')
    
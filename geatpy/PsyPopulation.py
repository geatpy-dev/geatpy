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
    
    Encodings : list  - 存储各染色体编码方式的列表。
    
    Fields    : list  - 存储各染色体对应的译码矩阵的列表。
    
    Chroms    : list  - 存储种群各染色体矩阵的列表。
    
    Linds     : list  - 存储种群各染色体长度的列表。
    
    ObjV      : array - 种群目标函数值矩阵，每一行对应一个个体的目标函数值，每一列对应一个目标。
    
    FitnV     : array - 种群个体适应度列向量，每个元素对应一个个体的适应度，最小适应度为0。
    
    CV        : array - CV(Constraint Violation Value)是用来定量描述违反约束条件程度的矩阵，每行对应一个个体，每列对应一个约束。
                        注意：当没有设置约束条件时，CV设置为None。
    
    Phen      : array - 种群表现型矩阵（即染色体解码后所代表的决策变量所组成的矩阵）。
    
函数:
    详见源码。

"""

    def __init__(self, Encodings, Fields, NIND, Chroms=None, ObjV=None, FitnV=None, CV=None, Phen=None):

        """
        描述: 种群类的构造函数，用于实例化种群对象，例如：
             import geatpy as ea
             population = ea.PsyPopulation(Encodings, Fields, NIND)，
             NIND为所需要的个体数，
             此时得到的population还没被真正初始化，仅仅是完成种群对象的实例化。
             该构造函数必须传入Chroms，才算是完成种群真正的初始化。
             一开始可以只传入Encodings, Fields以及NIND来完成种群对象的实例化，
             其他属性可以后面再通过计算进行赋值。
             
        """

        if isinstance(NIND, int) and NIND >= 0:
            self.sizes = NIND
        else:
            raise RuntimeError('error in PysPopulation: Size error. (种群规模设置有误，必须为非负整数。)')
        self.ChromNum = len(Encodings)
        if self.ChromNum == 1:
            raise RuntimeError(
                'error in PysPopulation: ChromNum must be bigger than 1. (使用PysPopulation类时，染色体数目必须大于1，否则应该使用Population类。)')
        self.Encodings = Encodings
        self.Fields = Fields.copy()
        self.Chroms = [None] * self.ChromNum  # 初始化Chroms为元素是None的列表
        self.Linds = []
        if Chroms is None:
            self.Linds = [0] * self.ChromNum
        else:
            for i in range(self.ChromNum):
                if Chroms[i] is not None:
                    self.Linds.append(Chroms[i].shape[1])
                    self.Chroms[i] = Chroms[i].copy() if Chroms[i] is not None else None
                else:
                    self.Linds.append(0)
        self.ObjV = ObjV.copy() if ObjV is not None else None
        self.FitnV = FitnV.copy() if FitnV is not None else None
        self.CV = CV.copy() if CV is not None else None
        self.Phen = Phen.copy() if Phen is not None else None

    def initChrom(self, NIND=None):

        """
        描述: 初始化种群染色体矩阵。

        输入参数:
            NIND : int - (可选参数)用于修改种群规模。
                         当其不缺省时，种群在初始化染色体矩阵前会把种群规模调整为NIND。

        输出参数:
            无输出参数。

        """

        if NIND is not None:
            self.sizes = NIND  # 重新设置种群规模
        # 遍历各染色体矩阵进行初始化
        for i in range(self.ChromNum):
            self.Chroms[i] = ea.crtpc(self.Encodings[i], self.sizes, self.Fields[i])  # 生成染色体矩阵
            self.Linds.append(self.Chroms[i].shape[1])  # 计算染色体的长度
        self.ObjV = None
        self.FitnV = np.ones((self.sizes, 1))  # 默认适应度全为1
        self.CV = None

    def decoding(self):

        """
        描述: 种群染色体解码。
        
        """

        Phen = np.ones((self.sizes, 0))  # 初始化一个空的矩阵
        # 遍历各染色体矩阵进行解码
        for i in range(self.ChromNum):
            if self.Encodings[i] == 'BG':  # 此时Field实际上为FieldD
                tempPhen = ea.bs2ri(self.Chroms[i], self.Fields[i])  # 把二进制/格雷码转化为实整数
            elif self.Encodings[i] == 'RI' or self.Encodings[i] == 'P':
                tempPhen = self.Chroms[i].copy()
            else:
                raise RuntimeError(
                    'error in PsyPopulation.decoding: Encoding must be ''BG'' or ''RI'' or ''P''. (编码设置有误，Encoding必须为''BG'', ''RI'' 或 ''P''。)')
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
        注意: index必须为一个slice或者为一个Numpy array类型的行向量或者为一个list类型的列表或者为一个整数，
             该函数不对传入的index参数的合法性进行更详细的检查。
        
        """

        # 计算切片后的长度
        if not isinstance(index, (slice, np.ndarray, list, int, np.int32, np.int64)):
            raise RuntimeError(
                'error in Population: index must be an integer, a 1-D list, or a 1-D array. (index必须是一个整数，一维的列表或者一维的向量。)')
        if isinstance(index, slice):
            NIND = (index.stop - (index.start if index.start is not None else 0)) // (
                index.step if index.step is not None else 1)
            index_array = index
        else:
            index_array = np.array(index).reshape(-1)
            if index_array.dtype == bool:
                NIND = int(np.sum(index_array))
            else:
                NIND = len(index_array)
            if len(index_array) == 0:
                index_array = []
        NewChroms = []
        for i in range(self.ChromNum):
            if self.Chroms[i] is None:
                raise RuntimeError('error in PsyPopulation: Chrom[i] is None. (种群染色体矩阵未初始化。)')
            NewChroms.append(self.Chroms[i][index_array])
        NIND = NewChroms[0].shape[0]
        return PsyPopulation(self.Encodings,
                             self.Fields,
                             NIND,
                             NewChroms,
                             self.ObjV[index_array] if self.ObjV is not None else None,
                             self.FitnV[index_array] if self.FitnV is not None else None,
                             self.CV[index_array] if self.CV is not None else None,
                             self.Phen[index_array] if self.Phen is not None else None)

    def shuffle(self):

        """
        shuffle : function - 打乱种群个体的个体顺序
        用法: 假设pop是一个种群矩阵，那么，pop.shuffle()即可完成对pop种群个体顺序的打乱。
        
        """

        shuff = np.arange(self.sizes)
        np.random.shuffle(shuff)  # 打乱顺序
        for i in range(self.ChromNum):
            if self.Chroms[i] is None:
                raise RuntimeError('error in PsyPopulation: Chrom[i] is None. (种群染色体矩阵未初始化。)')
            self.Chroms[i] = self.Chroms[i][shuff, :]
        self.ObjV = self.ObjV[shuff, :] if self.ObjV is not None else None
        self.FitnV = self.FitnV[shuff] if self.FitnV is not None else None
        self.CV = self.CV[shuff, :] if self.CV is not None else None
        self.Phen = self.Phen[shuff, :] if self.Phen is not None else None

    def __setitem__(self, index, pop):  # 种群个体赋值（种群个体替换）

        """
        描述: 种群个体的赋值
        用法: 假设pop是一个包含多于2个个体的种群矩阵，pop1是另一个包含2个个体的种群矩阵，那么
             pop[[0,1]] = pop1，即可完成将pop种群的第1、2个个体赋值为pop1种群的个体。
        注意: index必须为一个slice或者为一个Numpy array类型的行向量或者为一个list类型的列表或者为一个整数，
             该函数不对传入的index参数的合法性进行更详细的检查。
             此外，进行种群个体替换后，该函数不会对适应度进行主动重置，
             如果因个体替换而需要重新对所有个体的适应度进行评价，则需要手写代码更新种群的适应度。
             
        """

        # 对index进行格式处理
        if not isinstance(index, (slice, np.ndarray, list, int, np.int32, np.int64)):
            raise RuntimeError(
                'error in Population: index must be an integer, a 1-D list, or a 1-D array. (index必须是一个整数，一维的列表或者一维的向量。)')
        if isinstance(index, slice):
            index_array = index
        else:
            index_array = np.array(index).reshape(-1)
            if len(index_array) == 0:
                index_array = []
        for i in range(self.ChromNum):
            if self.Encodings[i] != pop.Encodings[i]:
                raise RuntimeError('error in PsyPopulation: Encoding disagree. (两种群染色体的编码方式必须一致。)')
            if np.all(self.Fields[i] == pop.Fields[i]) == False:
                raise RuntimeError('error in PsyPopulation: Field disagree. (两者的译码矩阵必须一致。)')
            if self.Chroms[i] is None:
                raise RuntimeError('error in PsyPopulation: Chrom[i] is None. (种群染色体矩阵未初始化。)')
            self.Chroms[i][index_array] = pop.Chroms[i]
        if self.ObjV is not None:
            if pop.ObjV is None:
                raise RuntimeError('error in PsyPopulation: ObjV disagree. (两者的目标函数值矩阵必须要么同时为None要么同时不为None。)')
            self.ObjV[index_array] = pop.ObjV
        if self.FitnV is not None:
            if pop.FitnV is None:
                raise RuntimeError('error in PsyPopulation: FitnV disagree. (两者的适应度列向量必须要么同时为None要么同时不为None。)')
            self.FitnV[index_array] = pop.FitnV
        if self.CV is not None:
            if pop.CV is None:
                raise RuntimeError('error in PsyPopulation: CV disagree. (两者的违反约束程度矩阵必须要么同时为None要么同时不为None。)')
            self.CV[index_array] = pop.CV
        if self.Phen is not None:
            if pop.Phen is None:
                raise RuntimeError('error in PsyPopulation: Phen disagree. (两者的表现型矩阵必须要么同时为None要么同时不为None。)')
            self.Phen[index_array] = pop.Phen
        self.sizes = self.Phen.shape[0]  # 更新种群规模

    def __add__(self, pop):

        """
        描述: 种群个体合并
        用法: 假设pop1, pop2是两个种群，它们的个体数可以相等也可以不相等，此时
             pop = pop1 + pop2，即可完成对pop1和pop2两个种群个体的合并。
        注意：
            进行种群合并后，该函数不会对适应度进行主动重置，
            如果因种群合并而需要重新对所有个体的适应度进行评价，则需要手写代码更新种群的适应度。
            
        """

        NIND = self.sizes + pop.sizes  # 得到合并种群的个体数
        NewChroms = self.Chroms
        for i in range(self.ChromNum):
            if self.Encodings[i] != pop.Encodings[i]:
                raise RuntimeError('error in PsyPopulation: Encoding disagree. (两种群染色体的编码方式必须一致。)')
            if np.all(self.Fields[i] == pop.Fields[i]) == False:
                raise RuntimeError('error in PsyPopulation: Field disagree. (两者的译码矩阵必须一致。)')
            if self.Chroms[i] is None or pop.Chroms[i] is None:
                raise RuntimeError('error in PsyPopulation: Chrom is None. (种群染色体矩阵未初始化。)')
            NewChroms[i] = np.vstack([NewChroms[i], pop.Chroms[i]])
        return PsyPopulation(self.Encodings,
                             self.Fields,
                             NIND,
                             NewChroms,
                             np.vstack(
                                 [self.ObjV, pop.ObjV]) if self.ObjV is not None and pop.ObjV is not None else None,
                             np.vstack(
                                 [self.FitnV, pop.FitnV]) if self.FitnV is not None and pop.FitnV is not None else None,
                             np.vstack([self.CV, pop.CV]) if self.CV is not None and pop.CV is not None else None,
                             np.vstack(
                                 [self.Phen, pop.Phen]) if self.Phen is not None and pop.Phen is not None else None)

    def __len__(self):

        """
        描述: 计算种群规模
        用法: 假设pop是一个种群，那么len(pop)即可得到该种群的个体数。
             实际上，种群规模也可以通过pop.sizes得到。
             
        """

        return self.sizes

    def save(self, dirName='Result'):

        """
        描述: 该函数将在字符串dirName所指向的文件夹下保存种群的信息，其中：
        "Encodingsi.txt"保存种群的染色体编码，i为0,1,2,3...；
        "Fieldsi.csv"保存种群染色体的译码矩阵，i为0,1,2,3...；
        "Chromsi.csv"保存种群的染色体矩阵，i为0,1,2,3...；
        "ObjV.csv"保存种群的目标函数矩阵；
        "FitnV.csv"保存种群个体的适应度列向量；
        "CV.csv"保存种群个体的违反约束程度矩阵；
        "Phen.csv"保存种群染色体表现型矩阵；
        注意：该函数不会对种群的合法性进行检查。
        
        """

        if self.sizes > 0:
            if os.path.exists(dirName) == False:
                os.makedirs(dirName)
            for i in range(self.ChromNum):
                with open(dirName + '/Encodings' + str(i) + '.txt', 'w') as file:
                    file.write(str(self.Encodings[i]))
                    file.close()
                np.savetxt(dirName + '/Fields' + str(i) + '.csv', self.Fields[i], delimiter=',')
                np.savetxt(dirName + '/Chroms' + str(i) + '.csv', self.Chroms[i], delimiter=',')
            if self.ObjV is not None:
                np.savetxt(dirName + '/ObjV.csv', self.ObjV, delimiter=',')
            if self.FitnV is not None:
                np.savetxt(dirName + '/FitnV.csv', self.FitnV, delimiter=',')
            if self.CV is not None:
                np.savetxt(dirName + '/CV.csv', self.CV, delimiter=',')
            if self.Phen is not None:
                np.savetxt(dirName + '/Phen.csv', self.Phen, delimiter=',')
            print('种群信息导出完毕。')

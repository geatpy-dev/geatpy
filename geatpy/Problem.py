# -*- coding: utf-8 -*-
import os
import numpy as np


class Problem:

    """
Problem : Class - 问题类

描述:
    问题类是用来存储与待求解问题相关信息的一个类。

属性:
    name      : str   - 问题名称（可以自由设置名称）。
    
    M         : int   - 目标维数，即有多少个优化目标。
    
    maxormins : array - 目标函数最小最大化标记的Numpy ndarray一维数组，1表示最小化，-1表示最大化，例如：
                        array([1,1,-1,-1])，表示前2个目标是最小化，后2个目标是最大化。
    
    Dim       : int   - 决策变量维数，即有多少个决策变量。
    
    varTypes  : array - 连续或离散标记，是Numpy ndarray类型的一维数组，
                        0表示对应的决策变量是连续的；1表示对应的变量是离散的。

    lb        : array - 存储着各个变量的下界。

    ub        : array - 存储着各个变量的上界。

    ranges    : array - 决策变量范围矩阵，第一行对应决策变量的下界，第二行对应决策变量的上界。
    
    borders   : array - 决策变量范围的边界矩阵，第一行对应决策变量的下边界，第二行对应决策变量的上边界，
                        0表示范围中不含边界，1表示范围包含边界。
    
    ReferObjV : array - 存储着目标函数参考值的矩阵，每一行对应一组目标函数参考值，每一列对应一个目标函数。

    TinyReferObjV : array - 从ReferObjV中均匀抽取的数目更少的目标函数参考值矩阵。

函数:
    aimFunc(pop)      : 目标函数，需要在继承类即自定义的问题类中实现，或是传入已实现的函数。
                        其中pop为Population类的对象，代表一个种群，
                        pop对象的Phen属性（即种群染色体的表现型）等价于种群所有个体的决策变量组成的矩阵，
                        该函数根据该Phen计算得到种群所有个体的目标函数值组成的矩阵，并将其赋值给pop对象的ObjV属性。
                        若有约束条件，则在计算违反约束程度矩阵CV后赋值给pop对象的CV属性（详见Geatpy数据结构）。
                        该函数不返回任何的返回值，求得的目标函数值保存在种群对象的ObjV属性中。
                        例如：population为一个种群对象，则调用aimFunc(population)即可完成目标函数值的计算，
                        此时可通过population.ObjV得到求得的目标函数值，population.CV得到违反约束程度矩阵。
                        注意：在子类中，aimFunc()和evalVars()两者只能重写一个。

    evalVars(v)       : 用于直接传入决策变量矩阵来计算对应的目标函数矩阵和违反约束程度矩阵。该函数需要被子类重写。

    evaluation(pop)   : 调用aimFunc()或evalVars()计算传入种群的目标函数值和违反约束程度。

    calReferObjV()    : 计算目标函数参考值，需要在继承类中实现，或是传入已实现的函数。
    
    getReferObjV()    : 获取目标函数参考值。

"""

    def __init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin=None, ubin=None, aimFunc=None, evalVars=None, calReferObjV=None):

        """
        描述:
            构造函数。
        """

        self.name = name
        self.M = M
        self.maxormins = None if np.all(maxormins == 1) or maxormins is None else np.array(maxormins)
        self.Dim = Dim
        self.varTypes = np.array(varTypes)
        self.lb = np.array(lb)
        self.ub = np.array(ub)
        self.ranges = np.array([lb, ub])  # 初始化ranges（决策变量范围矩阵）
        if lbin is None:
            lbin = [1] * Dim
        if ubin is None:
            ubin = [1] * Dim
        self.borders = np.array([lbin, ubin])  # 初始化borders（决策变量范围边界矩阵）
        self.aimFunc = aimFunc if aimFunc is not None else self.aimFunc  # 初始化目标函数接口
        self.evalVars = evalVars if evalVars is not None else self.evalVars
        self.calReferObjV = calReferObjV if calReferObjV is not None else self.calReferObjV  # 初始化理论最优值计算函数接口
        self.ReferObjV = self.getReferObjV()  # 计算目标函数参考值
        if self.ReferObjV is not None:
            if self.ReferObjV.shape[0] > 100:
                chooseIdx = np.linspace(0, self.ReferObjV.shape[0] - 1, 100).astype(np.int32)
                self.TinyReferObjV = self.ReferObjV[chooseIdx, :]
            else:
                self.TinyReferObjV = self.ReferObjV
        else:
            self.TinyReferObjV = None

    def aimFunc(self, pop):

        """
        描述:
            该函数需要被子类重写，用于计算整个种群的目标函数矩阵和违反约束程度矩阵。

        注意:
            如果执行到这里抛出了异常，说明自定义的问题类既没有重写aimFunc()，也没有重写evalVars()。
            aimFunc和evalVars两者中必须有一个被子类重写。
            若子类同时重写了aimFunc和evalVars，则evalVars无效。
        
        输入参数:
            pop : class <Population> - 种群对象

        输出参数:
            无输出参数。
        
        """

        raise RuntimeError('error in Problem: aim function has not been initialized. (未在问题子类中设置目标函数！)')

    def evalVars(self, Vars):

        """
        描述:
            evalVars即evaluate variables，用于直接传入决策变量矩阵来计算对应的目标函数矩阵和违反约束程度矩阵。该函数需要被子类重写。

        用法:
            ObjV = evalVars(Vars)
            ObjV, CV = evalVars(Vars)

        注意:
            如果执行到这里抛出了异常，说明自定义的问题类既没有重写aimFunc()，也没有重写evalVars()。
            aimFunc和evalVars两者的其中一个必须有一个被子类重写。
            若子类同时重写了aimFunc和evalVars，则evalVars无效。

        输入参数:
            Vars : array - 决策变量矩阵。每一行代表一组决策变量。

        输出参数:
            ObjV : array - 目标函数值矩阵。

            CV : array - 违反约束程度矩阵。

        """

        raise RuntimeError('error in Problem: aim function has not been initialized. (未在问题子类中设置目标函数！)')

    @staticmethod
    def single(func):
        """
        描述: 装饰器single。通过给目标函数添加装饰器，可以更专注于问题的模型本身。因为此时传入自定义目标函数的只有一个个体或一组决策变量。
        用法: 1. 给aimFunc()加上single装饰器标记后，可以让aimFunc()的传入参数：种群对象pop只有一个个体。
             2. 给evalVars()加上single装饰器标记后，可以让evalVars()的传入参数：Vars变成Numpy ndarray一维数组。即只传入一组决策变量。
        """
        def wrapper(param):
            if type(param) == np.ndarray:
                if param.ndim == 1:
                    param = np.atleast_2d(param)
                elif param.ndim > 2:
                    raise RuntimeError('Invalid input parameter of evalVars. (evalVars的传入参数非法，必须为Numpy ndarray一维数组或二维数组。)')
                return_object = func(param[0, :])
                if type(return_object) != tuple:
                    ObjV = [return_object]
                    for i in range(1, param.shape[0]):
                        Vars = param[i, :]
                        ObjV_i = func(Vars)
                        ObjV.append(np.atleast_1d(ObjV_i))
                    return np.vstack(ObjV)
                else:
                    ObjV_i, CV_i = return_object
                    ObjV = [ObjV_i]
                    CV = [CV_i]
                    for i in range(1, param.shape[0]):
                        Vars = param[i, :]
                        return_object = func(Vars)
                        ObjV_i, CV_i = return_object
                        ObjV.append(np.atleast_1d(ObjV_i))
                        CV.append(np.atleast_1d(CV_i))
                    return np.vstack(ObjV), np.vstack(CV)
            else:
                pop = param
                ObjV = []
                CV = []
                for i in range(param.sizes):
                    pop_i = param[i]
                    func(pop_i)
                    ObjV.append(pop_i.ObjV[0, :])
                    CV.append(pop_i.CV[0, :])
                pop.ObjV = np.vstack(ObjV)
                pop.CV = np.vstack(CV)

        return wrapper

    def evaluation(self, pop):

        """
        描述:
            调用aimFunc()或evalVars()计算传入种群的目标函数值和违反约束程度。
        注意:
            若子类同时重写了aimFunc和evalVars，则evalVars无效。

        输入参数:
            pop : class <Population> - 种群对象。

        """

        aimFuncCallFlag = False
        evalVarsCallFlag = False
        # 先尝试执行aimFunc()
        if self.aimFunc.__name__ != 'wrapper':  # 检查目标函数是否套了装饰器
            if self.aimFunc.__module__ != 'geatpy.Problem':
                aimFuncCallFlag = True
        else:  # 若套了装饰器，则aimFunc一定被重写了。
            aimFuncCallFlag = True
        if aimFuncCallFlag:
            self.aimFunc(pop)
        else:  # 再尝试执行evalVars()
            if self.evalVars.__name__ != 'wrapper':
                if self.evalVars.__module__ != 'geatpy.Problem':
                    evalVarsCallFlag = True
            else:  # 若套了装饰器，则evalVars一定被重写了。
                evalVarsCallFlag = True
        if evalVarsCallFlag:
            return_object = self.evalVars(pop.Phen)
            if type(return_object) != tuple:
                pop.ObjV = return_object
            else:
                pop.ObjV, pop.CV = return_object
        if not aimFuncCallFlag and not evalVarsCallFlag:
            raise RuntimeError('error in Problem: one of the function aimFunc and evalVars should be rewritten. (aimFunc和evalVars两个函数必须至少有一个被子类重写。)')

    def calReferObjV(self):

        """
        描述:
            如果待优化的模型知道理论全局最优解，则可以在自定义问题类里重写该函数，求出理论全局最优解对应的目标函数值矩阵。
            
        """

        return None

    def getReferObjV(self, reCalculate=False):

        """
        描述: 
            该函数用于读取/计算问题的目标函数参考值，这个参考值可以是理论上的全局最优解的目标函数值，也可以是人为设定的非最优的目标函数参考值。
            在获取或计算出理论全局最优解后，
            结果将被按照“问题名称_目标维数_决策变量个数.csv”的文件命名保存到referenceObjV文件夹内。
        
        输入参数:
            reCalculate : bool - 表示是否要调用calReferObjV()来重新计算目标函数参考值。
                                 当缺省时默认为False。

        输出参数:
            referenceObjV : array - 存储着目标函数参考值的矩阵，每一行对应一组目标函数参考值，每一列对应一个目标函数。
        
        """

        if not reCalculate:
            if self.calReferObjV.__module__ != 'geatpy.Problem':
                # 尝试读取数据
                if os.path.exists('referenceObjV'):
                    if os.path.exists('referenceObjV/' + self.name + '_M' + str(self.M) + '_D' + str(self.Dim) + '.csv'):
                        return np.atleast_2d(np.loadtxt('referenceObjV/' + self.name + '_M' + str(self.M) + '_D' + str(self.Dim) + '.csv',
                                          delimiter=','))
        # 若找不到数据，则调用calReferObjV()计算目标函数参考值
        referenceObjV = self.calReferObjV()
        if referenceObjV is not None:
            # 简单检查referenceObjV的合法性
            if not isinstance(referenceObjV, np.ndarray) or referenceObjV.ndim != 2 or referenceObjV.shape[1] != self.M:
                raise RuntimeError(
                    'error: ReferenceObjV is illegal. (目标函数参考值矩阵的数据格式不合法，请检查自定义问题类中的calReferObjV('
                    ')函数的代码，假如没有目标函数参考值，则在问题类中不需要定义calReferObjV()函数。)')
            # 保存数据
            if not os.path.exists('referenceObjV'):
                os.makedirs('referenceObjV')
            np.savetxt('referenceObjV/' + self.name + '_M' + str(self.M) + '_D' + str(self.Dim) + '.csv', referenceObjV,
                       delimiter=',')
        self.ReferObjV = referenceObjV
        return referenceObjV

    def __str__(self):
        info = {}
        info['name'] = self.name
        info['M'] = self.M
        info['maxormins'] = self.maxormins
        info['Dim'] = self.Dim
        info['varTypes'] = self.varTypes
        info['lb'] = self.lb
        info['ub'] = self.ub
        info['borders'] = self.borders
        return str(info)
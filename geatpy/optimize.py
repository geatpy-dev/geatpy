# -*- coding: utf-8 -*-

import os
import time
import numpy as np
import geatpy as ea


def optimize(algorithm,
             seed=None,
             prophet=None,
             verbose=None,
             drawing=None,
             outputMsg=True,
             drawLog=True,
             saveFlag=True,
             dirName=None,
             **kwargs):

    """

    描述: 该函数用于快速部署问题和算法，然后对问题进行求解。

    输入参数:
        algorithm : <class: class> - 算法类的引用。

        seed      : int  - 随机数种子。

        prophet   : <class: Population> / Numpy ndarray - 先验知识。可以是种群对象，
                                                          也可以是一组或多组决策变量组成的矩阵（矩阵的每一行对应一组决策变量）。

        verbose   : bool - 控制是否在输入输出流中打印输出日志信息。
                           该参数将被传递给algorithm.verbose。
                           如果algorithm已设置了该参数的值，则调用optimize函数时，可以不传入该参数。

        drawing   : int  - 算法类控制绘图方式的参数，
                           0表示不绘图；
                           1表示绘制最终结果图；
                           2表示实时绘制目标空间动态图；
                           3表示实时绘制决策空间动态图。
                           该参数将被传递给algorithm.drawing。
                           如果algorithm已设置了该参数的值，则调用optimize函数时，可以不传入该参数。

        outputMsg : bool - 控制是否输出结果以及相关指标信息。

        drawLog   : bool - 用于控制是否根据日志绘制迭代变化图像。

        saveFlag  : bool - 控制是否保存结果。

        dirName   : str  - 文件保存的路径。当缺省或为None时，默认保存在当前工作目录的'result of job xxxx-xx-xx xxh-xxm-xxs'文件夹下。

    输出参数:
        result    : dict - 一个保存着结果的字典。内容为：
                           {'success': True or False,  # 表示算法是否成功求解。
                            'stopMsg': xxx,  # 存储着算法停止原因的字符串。
                            'optPop': xxx,  # 存储着算法求解结果的种群对象。如果无可行解，则optPop.sizes=0。optPop.Phen为决策变量矩阵，optPop.ObjV为目标函数值矩阵。
                            'lastPop': xxx,  # 算法进化结束后的最后一代种群对象。
                            'Vars': xxx,  # 等于optPop.Phen，此处即最优解。若无可行解，则Vars=None。
                            'ObjV': xxx,  # 等于optPop.ObjV，此处即最优解对应的目标函数值。若无可行解，ObjV=None。
                            'CV': xxx,  # 等于optPop.CV，此处即最优解对应的违反约束程度矩阵。若无可行解，CV=None。
                            'startTime': xxx,  # 程序执行开始时间。
                            'endTime': xxx,  # 程序执行结束时间。
                            'executeTime': xxx,  # 算法所用时间。
                            'nfev': xxx,  # 算法评价次数
                            'gd': xxx,  # (多目标优化且给定了理论最优解时才有) GD指标值。
                            'igd': xxx,  # (多目标优化且给定了理论最优解时才有) IGD指标值。
                            'hv': xxx,  # (多目标优化才有) HV指标值。
                            'spacing': xxx}  # (多目标优化才有) Spacing指标值。

    """

    startTime = time.strftime("%Y-%m-%d %Hh-%Mm-%Ss", time.localtime())  # 记录程序开始时间
    if dirName is None:
        dirName = 'result of job ' + str(startTime)
    if dirName != '':
        dirName += '/'
    if saveFlag:
        if dirName != '':
            if not os.path.exists(dirName):
                os.makedirs(dirName)
        algorithm.dirName = dirName  # 只有在设置了saveFlag=True时，才让algorithm的dirName同步成一致。
    if seed is not None:
        np.random.seed(seed)
    # 处理先验知识
    prophetPop = None
    if prophet is not None:
        if type(prophet) == np.ndarray:
            if prophet.ndim != 2:
                prophet = np.atleast_2d(prophet)
            if algorithm.population.Encoding == 'RI' or algorithm.population.Encoding == 'P':
                prophetPop = ea.Population(algorithm.population.Encoding, algorithm.population.Field, 1, prophet)
            elif algorithm.population.Encoding == 'BG':
                prophetPop = ea.Population(algorithm.population.Encoding, algorithm.population.Field, 1, ea.ri2bs(prophet, algorithm.population.Field))
            else:
                raise RuntimeError('error in optimize: The encoding should be ''RI'', ''P'' or ''BG''. (种群编码必须为''RI'', ''P'' 或 ''BG''。)')
        elif type(prophet) == ea.Population:
            prophetPop = prophet.copy()
        else:
            raise RuntimeError('error in optimize: The type of prophet must be Numpy ndarray or Population. (prophet的类型必须为Numpy ndarray数组或种群类。)')
    # 参数设置
    algorithm.verbose = verbose if verbose is not None else algorithm.verbose
    algorithm.drawing = drawing if drawing is not None else algorithm.drawing
    # 开始求解
    [optPop, lastPop] = algorithm.run(prophetPop)
    # 生成结果
    result = {}
    result['success'] = True if optPop.sizes > 0 else False
    result['stopMsg'] = algorithm.stopMsg
    result['optPop'] = optPop
    result['lastPop'] = lastPop
    result['Vars'] = optPop.Phen if optPop.sizes > 0 else None
    result['ObjV'] = optPop.ObjV if optPop.sizes > 0 else None
    result['CV'] = optPop.CV if optPop.sizes > 0 else None
    result['executeTime'] = algorithm.passTime
    result['nfev'] = algorithm.evalsNum
    # 计算指标
    if algorithm.problem.M > 1 and optPop.sizes != 0:
        if algorithm.problem.ReferObjV is not None:
            GD = ea.indicator.GD(optPop.ObjV, algorithm.problem.ReferObjV)
            result['gd'] = GD
            IGD = ea.indicator.IGD(optPop.ObjV, algorithm.problem.ReferObjV)
            result['igd'] = IGD
        HV = ea.indicator.HV(optPop.ObjV)
        result['hv'] = HV
        Spacing = ea.indicator.Spacing(optPop.ObjV)
        result['spacing'] = Spacing
    # 输出结果
    if outputMsg:
        print('Execution time: %s s' % algorithm.passTime)
        print('Evaluation number: %s' % algorithm.evalsNum)
        if algorithm.problem.M == 1:
            if optPop.sizes != 0:
                print('The best objective value is: %s' % optPop.ObjV[0][0])
                print('The best variables are: ')
                for i in range(optPop.Phen.shape[1]):
                    print(optPop.Phen[0, i], end='\t')
                print()
            else:
                print('Optimization fail: Could not find any feasible solution.')
        else:
            if optPop.sizes != 0:
                print('The number of non-dominated solutions is: %d' % optPop.sizes)
            if optPop.sizes != 0:
                if algorithm.problem.ReferObjV is not None:
                    print('gd: %.5f' % result['gd'])
                    print('igd: %.5f' % result['igd'])
                print('hv: %.5f' % result['hv'])
                print('spacing: %.5f' % result['spacing'])
            else:
                print('Optimization fail: Could not find any feasible solution.')
    # 根据算法执行日志绘制迭代变化曲线
    if drawLog and algorithm.log is not None:
        if algorithm.problem.M == 1:
            trace = np.array(algorithm.log['f_opt'])
            plotter = ea.ParCoordPlotter(len(trace), grid=True, legend=True, title='Best Objective Value Trace Plot', coordLabels=['Generation Number', 'Value'], saveName=dirName + 'Best Objective Value Trace Plot' if saveFlag else None)
            plotter.add(trace, color='blue', label='Best-found objective value')
            plotter.draw()
            plotter.show()
        else:
            drawNameList = ['GD', 'IGD', 'HV', 'Spacing']
            for drawName in drawNameList:
                trace = np.array(algorithm.log[drawName.lower()])
                if trace[0] is not None:
                    plotter = ea.ParCoordPlotter(len(trace), grid=True, legend=True, title=drawName + ' Trace Plot', coordLabels=['Generation Number', 'Value'], saveName=dirName + drawName + ' Trace Plot' if saveFlag else None)
                    plotter.add(trace, color='blue', label=drawName)
                    plotter.draw()
                    plotter.show()
    endTime = time.strftime("%Y-%m-%d %Hh-%Mm-%Ss", time.localtime())  # 记录程序结束时间
    result['startTime'] = startTime
    result['endTime'] = endTime
    # 保存结果
    if saveFlag:
        optPop.save(dirName+'optPop')
        # 记录问题
        with open(dirName + 'problem info.txt', 'w') as file:
            file.write(str(algorithm.problem))
        # 记录种群参数
        popInfo = {}
        if isinstance(algorithm.population, ea.Population):
            popInfo['Population Info'] = algorithm.population.getInfo()
        elif isinstance(algorithm.population, list):
            infoList = []
            for pop in algorithm.population:
                infoList.append(pop.getInfo())
            popInfo['Population Info'] = infoList
        with open(dirName + 'population info.txt', 'w') as file:
            file.write(str(popInfo))
        # 记录算法参数
        with open(dirName + 'algorithm info.txt', 'w') as file:
            file.write(str(algorithm))
        # 记录结果
        resultInfo = result.copy()
        resultInfo.pop('optPop')
        resultInfo.pop('lastPop')
        with open(dirName + 'result info.txt', 'w') as file:
            file.write(str(resultInfo))

    return result

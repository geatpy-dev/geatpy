# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation


class PointScatter:
    """
    PointScatter: 绘制2维或3维散点图。

    属性:
        Dimension   : int - 数据点的维度。它必须等于add(points)中points列数(points为一个2D的Numpy ndarray数组)。

        grid        : bool - 控制是否绘制网格。

        legend      : bool - 控制是否绘制图例。

        title       : str  - 图片标题。

        coordLabels : list - 存储着各个坐标轴名称的列表。

        saveName    : str  - 保存图片的文件名（不含文件格式后缀。保存的文件为.svg格式）
                             当缺省或为None时，不保存图片。
    
    """

    def __init__(self, Dimension, grid=True, legend=True, title=None, coordLabels=None, saveName=None):
        if title is None:
            title = 'ParCoordPlot'
        if coordLabels is None:
            coordLabels = ['F1', 'F2', 'F3']
        # check Dimension
        if Dimension != 2 and Dimension != 3:
            raise RuntimeError('Error in PointScatter.py: M must be 2 or 3. (M 必须为2或3。)')
        self.Dimension = Dimension
        self.grid = grid
        self.legend = legend
        self.title = title
        self.coordLabels = coordLabels
        self.saveName = saveName
        self.data_set = []  # 存储用于绘图的数据，每次调用addPoints，都会往data_set里添加一组绘图数据。
        self.params_set = []  # 存储所有的绘图参数。其每个元素是一个dict，格式为：{'marker': xxx, 'objectsize': xxx, 'color': xxx, 'alpha': xxx, 'label': xxx}。
        self.plots_set = []  # 存储图形中的所有图形对象。
        self.history = []  # 用于在绘制动画时保存历史数据。
        self.fig = None
        self.ax = None

    def add(self, points, marker='o', objectSize=5, color='blue', alpha=1.0, label=None):
        """
        用于添加一组点到图形中。输入参数points的每一行表示每个点的坐标值。
        当输入参数points只有一行，或者是1D的Numpy ndarray一维数组时，只添加一个点。
        objectSize控制点的大小以及线的粗细。
        """
        if points is None:
            return
        if not isinstance(points, np.ndarray):
            raise RuntimeError('Error in PointScatter.py: The type of the points must be numpy ndarray. (points必须是Numpy ndarray类型。)')
        if points.ndim == 1:
            points = np.array([points])
            if points.shape[1] != self.Dimension:
                raise RuntimeError(
                    'Error in PointScatter.py: The length of the points must be equal to Dimension if its dimension is 1. (points是1维时，其长度必须等于Dimension。)')
        elif points.ndim == 2:
            if points.shape[0] == 0:
                return
            if points.shape[1] != self.Dimension:
                raise RuntimeError(
                    'Error in PointScatter.py: The number of the column of the points must be equal to Dimension if its dimension is 2. ('
                    'points是2维时，其列数必须等于Dimension。)')
        self.data_set.append(points)
        self.params_set.append({'marker': marker, 'objectsize': objectSize, 'color': color, 'alpha': alpha, 'label': label})
        self.history += [{'data_set': self.data_set, 'params_set': self.params_set}]

    def draw(self):
        # 开始绘制
        if self.Dimension == 2:
            if self.fig is None and self.ax is None:
                self.fig, self.ax = plt.subplots()  # 生成一块画布和创建绘图区域
            for idx, data in enumerate(self.data_set):
                params = self.params_set[idx]
                plot = self.ax.plot(data[:, 0], data[:, 1], params['marker'], markersize=params['objectsize'], linewidth=params['objectsize'], color=params['color'], alpha=params['alpha'], label=params['label'])
                self.plots_set.append(plot)
                self.ax.set_xlabel(self.coordLabels[0])
                self.ax.set_ylabel(self.coordLabels[1])
        elif self.Dimension == 3:
            if self.fig is None and self.ax is None:
                self.fig = plt.figure()  # 生成一块画布
                self.ax = self.fig.add_subplot(111, projection='3d')  # 创建绘图区域
                self.ax.view_init(elev=30, azim=45)  # 旋转
            for idx, data in enumerate(self.data_set):
                params = self.params_set[idx]
                plot = self.ax.plot(data[:, 0], data[:, 1], data[:, 2], params['marker'], markersize=params['objectsize'], color=params['color'], alpha=params['alpha'], label=params['label'])
                self.plots_set.append(plot)
                self.ax.set_xlabel(self.coordLabels[0])
                self.ax.set_ylabel(self.coordLabels[1])
                self.ax.set_zlabel(self.coordLabels[2])
        if self.title is not None:
            self.ax.set_title(self.title)
        if self.legend:
            plt.legend()
        plt.grid(self.grid)
        plt.draw()

    def refresh(self):
        if self.fig and self.ax:
            plt.pause(1/24)
            self.ax.cla()
        self.data_set = []
        self.params_set = []
        self.plots_set = []

    def show(self):
        if self.saveName is not None:
            self.fig.savefig(self.saveName + '.svg', dpi=300, bbox_inches='tight')
        plt.show()

    def createAnimation(self, fps=6):
        """
        该函数根据self.history记录的数据，绘制动画并保存到文件中。
        fps表示每秒钟绘制多少帧。
        """
        def update(i, plotObject):
            plotObject.ax.cla()
            plotObject.data_set = plotObject.history[i]['data_set']
            plotObject.params_set = plotObject.history[i]['params_set']
            plotObject.draw()
        if len(self.history) > 0:
            if self.fig is None and self.ax is None:
                if self.Dimension == 2:
                    self.fig, self.ax = plt.subplots()  # 生成一块画布和创建绘图区域
                elif self.Dimension == 3:
                    self.fig = plt.figure()  # 生成一块画布
                    self.ax = Axes3D(self.fig)  # 创建绘图区域
                    self.ax.view_init(elev=30, azim=45)  # 旋转
            print('Creating gif...')
            anim = animation.FuncAnimation(self.fig, update, frames=len(self.history), fargs=(self, ))
            anim.save(self.title + '.gif', fps=fps)
            print('gif has been created.')

    def close(self):
        plt.close()
        for item in self.history:
            item.clear()
        self.history = []
        self.data_set = []
        self.params_set = []
        self.plots_set = []

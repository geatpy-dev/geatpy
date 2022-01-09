# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class ParCoordPlotter:
    """
    ParCoordPlotter: 绘制平行坐标图。

    属性:
        Dimension   : int  - 数据点的维度。它必须等于add(points)中points列数(points为一个2D的Numpy ndarray数组)。

        xtickList   : list - 存储着平行坐标图横坐标数据的列表。

        grid        : bool - 控制是否绘制网格。

        legend      : bool - 控制是否绘制图例。

        title       : str  - 图片标题。

        coordLabels : list - 存储着各个坐标轴名称的列表。

        saveName    : str  - 保存图片的文件名（不含文件格式后缀。保存的文件为.svg格式）
                             当缺省或为None时，不保存图片。

    """

    def __init__(self, Dimension, xtickList=None, grid=True, legend=True, title=None, coordLabels=None, saveName=None):
        if title is None:
            title = 'ParCoordPlot'
        if coordLabels is None:
            coordLabels = ['Dimension Number', 'Value']
        self.Dimension = Dimension
        self.xtickList = xtickList
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

    def add(self, points, marker='-', objectSize=2, color='blue', alpha=1.0, label=None):
        """
        用于添加一组点到图形中。输入参数points的每一行表示每个点的坐标值。
        当输入参数points只有一行，或者是1D的Numpy ndarray一维数组时，只添加一个点。
        objectSize控制点的大小以及线的粗细。
        """

        if points is None:
            return
        if not isinstance(points, np.ndarray):
            raise RuntimeError('Error in ParCoordPlotter.py: The type of the points must be numpy ndarray. (points必须是Numpy ndarray类型。)')
        if points.ndim == 1:
            points = np.array([points])
            if points.shape[1] != self.Dimension:
                raise RuntimeError(
                    'Error in ParCoordPlotter.py: The length of the points must be equal to Dimension if its dimension is 1. (points是1维时，其长度必须等于Dimension。)')
        elif points.ndim == 2:
            if points.shape[0] == 0:
                return
            if points.shape[1] != self.Dimension:
                raise RuntimeError(
                    'Error in ParCoordPlotter.py: The number of the column of the points must be equal to Dimension if its dimension is 2. ('
                    'points是2维时，其列数必须等于Dimension。)')
        self.data_set.append(points)
        self.params_set.append({'marker': marker, 'objectsize': objectSize, 'color': color, 'alpha': alpha, 'label': label})
        self.history += [{'data_set': self.data_set, 'params_set': self.params_set}]

    def draw(self):
        # 开始绘制
        if self.fig is None and self.ax is None:
            self.fig, self.ax = plt.subplots()  # 生成一块画布和创建绘图区域
        for idx, data in enumerate(self.data_set):
            params = self.params_set[idx]
            x = np.tile(np.arange(1, data.shape[1] + 1).reshape(-1, 1), (1, data.shape[0]))
            _data = data.T
            self.ax.plot(x[:, 0], _data[:, 0], params['marker'], markersize=params['objectsize'], linewidth=params['objectsize'], color=params['color'], alpha=params['alpha'], label=params['label'])
            plot = self.ax.plot(x, _data, params['marker'], markersize=params['objectsize'], linewidth=params['objectsize'], color=params['color'], alpha=params['alpha'])
            self.plots_set.append(plot)
        self.ax.set_xlabel(self.coordLabels[0])
        self.ax.set_ylabel(self.coordLabels[1])
        if self.xtickList is not None:
            plt.xticks(np.arange(1, len(self.xtickList) + 1), self.xtickList)
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
                self.fig, self.ax = plt.subplots()  # 生成一块画布和创建绘图区域
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

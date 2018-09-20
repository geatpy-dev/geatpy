import numpy as np

"""目标函数"""

def aimfuc(Phen):
    x1 = Phen[:, 0]
    x2 = Phen[:, 1]
    f = 21.5 + x1 * np.sin(4 * np.pi * x1) + x2 * np.sin(20 * np.pi * x2)
    return np.array([f]).T

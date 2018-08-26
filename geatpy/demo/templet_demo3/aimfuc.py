import numpy as np

def aimfuc(Phen):
    x1 = Phen[:, 0]
    x2 = Phen[:, 1]
    x3 = Phen[:, 2]
    x4 = Phen[:, 3]
    f = 18 * x1 + 10 * x2 + 12 * x3 + 8 * x4
    # 约束条件
    f[np.where(12 * x1 + 6 * x2 + 10 * x3 + 4 * x4 > 20)[0]] = 0
    f[np.where(x3 + x4 > 1)[0]] = 0
    f[np.where(x3 - x1 > 0)[0]] = 0
    f[np.where(x4 - x2 > 0)[0]] = 0
    return np.array([f]).T

# -*- coding: utf-8 -*-
"""
罚函数示例

Created on Sun Jul 29 15:17:10 2018

@author: jazzbin
"""

import numpy as np

def punishing(Phen, FitnV):
    idx = np.where(Phen == 0)[0]
    FitnV[idx] = np.min(FitnV) // 2
    return FitnV

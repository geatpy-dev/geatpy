# -*- coding: utf-8 -*-
"""
execScript: main.py
description:
    Use GA algorithm to search the sentence: Tom is a little boy.
"""
import numpy as np
import geatpy as ga

strs = 'Tom is a little boy.' # define the sentence
words = []
for c in strs:
    words.append(ord(c)) # change the words to ascii code.

# Attention: you had better put the aim-function in another file.
def aim(Phen): # define the aim function
    real = words
    diff = np.sum((Phen - real)**2, 1)
    return np.array([diff]).T

if __name__ == "__main__":
    AIM_M = __import__('main') # get the handle of aim-function
    # variables setting
    ranges = np.vstack([32 * np.ones((1, len(words))), 122 * np.ones((1, len(words)))]) # ranges of variables
    borders = np.vstack([np.ones((1, len(words))), 122 * np.ones((1, len(words)))]) # borders of variables
    FieldDR = ga.crtfld(ranges, borders) # create FieldDR
    # call the GEA algorithm template
    [pop_trace, var_trace, times] = ga.sga_new_real_templet(AIM_M, 'aim', None, None, FieldDR, problem = 'I', maxormin = 1, MAXGEN = 2000, NIND = 50, SUBPOP = 1, GGAP = 0.9, selectStyle = 'sus', recombinStyle = 'xovdp', recopt = 0.9, pm = 0.1, drawing = 1)
    # output results
    for num in var_trace[np.argmin(pop_trace[:, 1]), :]:
        print(chr(int(num)), end = '')
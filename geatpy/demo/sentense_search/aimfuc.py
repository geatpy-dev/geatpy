import numpy as np

def aimfuc(Phen):
    real = np.array([ord('I'),ord(' '),ord('a'),ord('m'),ord(' '),ord('a'),
                     ord(' '),ord('l'),ord('i'),ord('t'),ord('t'),ord('l'),
                     ord('e'),ord(' '),ord('b'),ord('o'),ord('y')])
    diff = np.sum((Phen - real)**2, 1)
    return np.array([diff]).T

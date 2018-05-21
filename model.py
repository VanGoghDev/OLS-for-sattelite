import numpy as np
import math
from abc import *


def getsize(a):
    size = 0
    for i in range(len(a)):
        size += 1
    return size


# noinspection PyCompatibility
class Model(metaclass=ABCMeta):
    # __metaclass__ = ABCMeta

    def __init__(self, at0, at1, ah, n):
        self.t0 = at0
        self.t1 = at1
        self.n = n
        self.sampling_increment = ah
        self.x0 = np.empty(n)
        self.result = np.empty((0, 0))
        self.set_result()

    def set_result(self):
        size = int(math.floor(self.t1 - self.t0))
        self.result = np.empty((size, self.n + 1))

    def add_result(self, a, t):
        k = t / self.sampling_increment
        k = int(math.floor(k))
        self.result[k][0] = t
        for i in range(1, getsize(self.x0) + 1):
            self.result[k][i] = a[i-1]

    def get_order(self):
        return getsize(self.x0)

    @abstractmethod
    def get_right(self, tv, t):
        pass

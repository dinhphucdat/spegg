import numpy as np

def getInt1D():
    return np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=np.int32)

def getInt2D():
    return np.array([
    [0, 1, 2, 3, 4],
    [5, 6, 7, 8, 9],
    [10, 11, 12, 13, 14]
    ], dtype=np.int32)

def getFloat1D():
    return np.array([0.0, 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9], dtype=np.float32)

def getFloat2D():
    return np.array([
    [0.0, 1.1, 2.2, 3.3, 4.4],
    [5.5, 6.6, 7.7, 8.8, 9.9],
    [10.1, 11.11, 12.12, 13.13, 14.14]
    ], dtype=np.float32)

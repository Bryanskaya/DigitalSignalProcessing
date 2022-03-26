import numpy.random
from numpy import sqrt
from scipy import signal
import numpy as np
import math
from random import random
import matplotlib.pyplot as plt
from scipy.signal import filtfilt, butter
from scipy.optimize import root_scalar
from numpy.fft import fft, ifft


def rhofunc(x, u1: np.array, u2: np.array, step: float, T: float, d: float, e:float):
    m = np.array(range(0, len(u1)), dtype=float)
    mult = step / len(u1)
    squ = 1 + (2 * math.pi / T * m) ** 2

    beta = mult * np.sum(x ** 2 * squ * abs(u1) ** 2 / (abs(u2) ** 2 * step ** 2 + x * squ) ** 2)
    gamma = mult * np.sum(abs(u2) ** 2 * step ** 2 * abs(u1) ** 2 * squ /
                       (abs(u2) ** 2 * step ** 2 + x * (1 + 2 * math.pi * m / T) ** 2) ** 2)
    y = beta - (d + e * np.sqrt(gamma)) ** 2
    return y


def tikhonov_f(u1: np.array, u2: np.array, step: float, T: float, d: float, e:float):
    func = lambda x: rhofunc(x, u1, u2, step, T, d, e)
    alpha = root_scalar(func, bracket=[0, 1]).root

    m = np.array(range(0, len(u1)), dtype=float)
    mult = step / len(u1)
    squ = 1 + (2 * math.pi / T * m) ** 2

    h = np.array(range(0, len(u1)), dtype=float)
    for k in range(len(h)):
        temp = np.exp(2 * math.pi * 1j * k * m / len(u1)) * u1 * np.conj(u2) / (abs(u2)**2 * step**2 + alpha * squ)
        h[k] = mult * np.sum(temp)
    return h


def main():
    s1 = np.sqrt(2.0)
    s2 = 1.0
    mult = 6
    step = 0.01
    NS = 0.01

    t = []
    i = -mult
    while i < mult + 1e-5:
        t.append(i)
        i += step

    # Сигнал Гаусса
    u1 = np.array(signal.gaussian(len(t), s1 / step))
    u2 = np.array(signal.gaussian(len(t), s2 / step))

    n1 = np.random.normal(-NS, NS, size=len(t))
    n2 = np.random.normal(-NS, NS, size=len(t))

    delta = np.std(n1, ddof=1)
    epsilon = np.std(n2, ddof=1)

    # Гаусс + шум
    x1 = u1 + n1
    x2 = u2 + n2

    # БПФ
    v1 = fft(x1)
    v2 = fft(x2)

    # регуляризация Тихонова
    filt = list(map(abs, ifft(fft(x2) * tikhonov_f(v1, v2, step, 2 * mult, delta, epsilon))))

    plt.plot(t, x1, label='Сигнал 1 + шум')
    plt.plot(t, x2, label='Сигнал 2 + шум')
    plt.plot(t, filt, label='C фильтром Тихонова')
    plt.legend()

    plt.show()


if __name__ == '__main__':
    main()

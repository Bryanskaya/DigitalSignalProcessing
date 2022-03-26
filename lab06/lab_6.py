from scipy import signal
from scipy.signal import filtfilt, butter
import numpy as np
import math
from random import random
import matplotlib.pyplot as plt


def impulse_noise(size, N, mult):
    step = math.floor(size/N)
    y = [0] * size
    for i in range(1, math.floor(N/2)):
        y[round(size/2) + i * step] = mult * (0.5 + random())
        y[round(size/2) - i * step] = mult * (0.5 + random())

    return y


def wiener(x, n):
    return [(1 - (ni / xi) ** 2) for xi, ni in zip(x, n)]


def main():
    sigma = 0.5

    n = 500
    step = 0.01
    #  n = int(input('Введите количество точек: '))
    #  step = float(input('Введите шаг: '))

    t_max = step * (n - 1) / 2
    t = []
    i = -t_max
    while i < t_max + 1e-5:
        t.append(i)
        i += step

    # Сигнал Гаусса
    x0 = list(signal.gaussian(len(t), sigma / step))

    NA = 0
    NS = 0.05
    n1 = []
    for _ in range(len(x0)):
        n1.append(np.random.normal(NA, NS))

    # Гуассовы помехи
    x1 = []
    for i in range(len(x0)):
        x1.append(x0[i] + n1[i])

    # Импульсные помехи
    count = 20
    M = 0.3
    n2 = impulse_noise(len(x0), count, M)
    x2 = []
    for i in range(len(x0)):
        x2.append(x0[i] + n2[i])

    y1 = wiener(np.fft.fft(x1), np.fft.fft(n1))
    y2 = wiener(np.fft.fft(x2), np.fft.fft(n2))

    graphics = [[x0, x1, x2],
                [x0, np.fft.ifft(np.fft.fft(x1) * y1), np.fft.ifft(np.fft.fft(x2) * y2)]]
    title1 = ['Сигнал Гаусса', 'Гауссовы помехи', 'Импульсные помехи']
    fig, axs = plt.subplots(2, 3)

    for i in range(len(graphics)):
        for j in range(len(graphics[0])):
            axs[i][j].plot(t, list(graphics[i][j]))
            axs[i][j].grid(True)
            axs[0][j].set_title(title1[j])
    axs[1][1].set_title('После фильтра Винера')

    plt.show()

if __name__ == '__main__':
    main()

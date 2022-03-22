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


def gaussian_filter_high(sigma, size):
    x = np.linspace(-size / 2, size / 2, size)
    y = []
    for i in range(len(x)):
        temp = -x[i] ** 2 / (2 * sigma ** 2)
        y.append(1 - math.exp(temp))

    s = sum(y)
    for i in range(len(y)):
        y[i] /= s

    return y


def butterworth_filter_high(D, size):
    x = np.linspace(-size / 2, size / 2, size)
    y = []
    for i in range(len(x)):
        y.append(1 / (1 + (D / x[i]) ** 4))
    s = sum(y)
    for i in range(len(y)):
        y[i] /= s

    return y


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
    count = 15
    M = 0.4
    n2 = impulse_noise(len(x0), count, M)
    x2 = []
    for i in range(len(x0)):
        x2.append(x0[i] + n2[i])

    G = gaussian_filter_high(4, 20)
    B = butterworth_filter_high(4, 20)
    #b, a = butter(6, 0.5)

    graphics = [x0, x1, x2]
    title1 = ['Сигнал Гаусса', 'Гауссовы помехи', 'Импульсные помехи']
    fig, axs = plt.subplots(3, 3)

    for i in range(len(graphics)):
        axs[0][i].plot(t, graphics[i])
        axs[0][i].grid(True)
        axs[0][i].set_title(title1[i])
    for i in range(len(graphics)):
        axs[1][i].plot(t, graphics[i] - filtfilt(G, 1, graphics[i]))
        axs[1][i].grid(True)
        axs[1][1].set_title('Фильтр низкий частот: Гаусса')
    for i in range(len(graphics)):
        axs[2][i].plot(t, graphics[i] - filtfilt(B, 1, graphics[i]))
        axs[2][i].grid(True)
        axs[2][1].set_title('Фильтр низкий частот: Баттерворт')

    plt.show()

if __name__ == '__main__':
    main()

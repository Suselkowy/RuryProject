import math
import matplotlib.pyplot as plt

import scipy.constants
import numpy as np

W1 = 1
X1 = 1/math.sqrt(3)
X2 = -X1

def change_to_curr_range(x, start_x, end_x):
    return (end_x - start_x) / 2 * x + (start_x + end_x) / 2


def base_function(i, h):
    def inside_base(x):
        if ((i - 1) * h < x <= i * h):
            return (x - h * (i - 1)) / h
        if (i * h < x < (i + 1) * h):
            return (h * (i + 1) - x) / h
        return 0

    return inside_base


def base_function_derivative(i, h):
    def inside_derevative(x):
        if ((i - 1) * h < x <= i * h):
            return 1 / h
        if (i * h < x < (i + 1) * h):
            return -1 / h
        return 0

    return inside_derevative


def combine_functions(f1, f2, negative=1):
    def kombinat(x):
        return f1(x) * f2(x) * negative

    return kombinat


def integral(function, start_x, end_x):
    x1 = change_to_curr_range(X1, start_x, end_x)
    x2 = change_to_curr_range(X2, start_x, end_x)
    w1 = W1

    return (end_x-start_x)/2 * (w1 * function(x1) + w1 * function(x2))

def p(x):
    return 10**11 if 1 < x <= 2 else 0

def L_p1(i, h):
    def inside(x):
        return p(x) * base_function(i, h)(x)
    return inside

def make_function(X, n, h):
    def inside(x):
        res = -5 * base_function(0, h)(x) - 4 * base_function(n, h)(x)
        for i in range(1, n):
            res += X[i-1] * base_function(i , h)(x)
        return res
    return inside

def plot(x,y):
    plt.plot(x,y)
    plt.show()

def start(n):
    matrix_A = np.zeros((n-1,n-1))
    matrix_B = np.zeros((n-1,1))
    h = 3 / n
    f1 = combine_functions(base_function_derivative(1, h), base_function_derivative(1, h), -1)
    matrix_A[0][0] = integral(f1, 0 * h, 1 * h) + integral(f1, 1 * h, 2 * h)

    for i in range(1,n):

        matrix_A[i-1][i-1] = matrix_A[0][0]

        if(i != n-1):
            f2 = combine_functions(base_function_derivative(i, h), base_function_derivative(i+1, h), -1)
            matrix_A[i-1][i] = integral(f2, i * h, (i+1) * h)
            matrix_A[i][i-1] = matrix_A[i-1][i]

        matrix_B[i-1][0] = 4 * math.pi * scipy.constants.G * (integral(L_p1(i, h), (i - 1) * h, i * h) + integral(L_p1(i, h), i * h, (i + 1) * h))

        l1 = combine_functions(base_function_derivative(0,h), base_function_derivative(i,h), -1)
        l2 = combine_functions(base_function_derivative(i,h), base_function_derivative(n,h), -1)
        matrix_B[i-1][0] += 5*integral(l1, (i-1) * h, i * h) + 4*integral(l2, i * h, (i+1) * h)

    x = np.linalg.solve(matrix_A, matrix_B)

    func = make_function(x, n, h)
    a = [h * i for i in range(n + 1)]
    b = [func(x) for x in a]
    plot(a,b)

def interface():
    try:
        n = int(input("Insert N: "))
        start(n)
    except ValueError:
        print("Provided value should be integer")
    except:
        print("unknown error occured")


interface()


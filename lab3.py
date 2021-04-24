from polynomial_algos import *
from matrix_algos import *
import system_algos

def task1():
    print('Первая ', end='')
    a = read_matrix()
    print('Вторая ', end='')
    b = read_matrix()
    print('Произведение:')
    print_matrix(mat_mul(a, b))

print('1. Умножение квадратных матриц разложением на блоки')
print('Пункт: ', end='')
[task1][int(input())-1]()

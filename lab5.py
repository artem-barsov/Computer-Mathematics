from polynomial_algos import *
from matrix_algos import *

def task1():
    print('Матрица:')
    a = read_matrix()
    print(polynom2str(char_polyn_LeVerrier(a), 'λ'))

def task2():
    print('Матрица:')
    a = read_matrix()
    print('Максимальное по модулю собственное значение λ =', eigenvalue(a))

print('1. Нахождение характеристического полинома произвольной квадратной матрицы по методу Леверрье')
print('2. Нахождение собственных значений матрицы')
print('Пункт: ', end='')
[task1, task2][int(input())-1]()

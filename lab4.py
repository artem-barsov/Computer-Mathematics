from polynomial_algos import *
from matrix_algos import *
import system_algos

def task1():
    print('Первая матрица:')
    a = read_matrix()
    print('Вторая матрица:')
    b = read_matrix()
    print('Произведение:')
    print_matrix(mat_mul(a, b))

def task2():
    print('Матрица:')
    a = read_matrix()
    try: a_inv = inverse_matrix(a)
    except Exception as e: print(e)
    else:
        print('Обратная матрица:')
        print_matrix(a_inv)
        print('Проверка умножением A * A^(-1):')
        print_matrix(mat_mul(a, a_inv))

def task3():
    print('Матрица:')
    a = read_matrix()
    try:
        print('Определитель матрицы =', '%g'%round(det(a), 3))
    except Exception as e:
        print(e)

print('1. Умножение квадратных матриц разложением на блоки')
print('2. Нахождение обартной матрицы через окаймляющие блоки')
print('3. Нахождение определителя матрицы через элементарные преобразования')
print('Пункт: ', end='')
[task1, task2, task3][int(input())-1]()

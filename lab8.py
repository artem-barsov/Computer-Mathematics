from system_algos import *
from polynomial_algos import get_script

def task1():
    a = [ [ 2, 7, -8, 6 ],
          [ 4, 4, 0, -7 ],
          [ -1,-3, 6, 3 ],
          [ 9,-7, -2,-8 ] ]
    b = [ -39, 41, 4, 113 ]
    print_system(a, b)
    x = Gauss(a, b)
    for i in range(1, len(x)+1):
        print(f'x{get_script(str(i), "sub")} = %g' % round(x[i-1], 3))

def task2():
    a = [ [ 10, 5,  0, 0, 0],
          [ 3, 10, -2, 0, 0],
          [ 0, 2, -9, -5, 0],
          [ 0, 0, 5, 16, -4],
          [ 0, 0, 0, -8, 16] ]
    b = [ -120, -91, 5, -74, -56 ]
    print_system(a, b)
    x = tridiagonal_algo(a, b)
    for i in range(1, len(x)+1):
        print(f'x{get_script(str(i), "sub")} = %g' % round(x[i-1], 3))

def task3():
    a = [ [ 24,   2,  4,  -9 ],
          [ -6, -27, -8,  -6 ],
          [ -4,   8, 19,   6 ],
          [  4,  -5, -3, -13 ] ]
    b = [ -9, -76, -79, -70 ]
    print_system(a, b)
    x = simple_iteration(a, b)
    for i in range(1, len(x)+1):
        print(f'x{get_script(str(i), "sub")} = %g' % round(x[i-1], 3))

def task4():
    a = [ [ 24,   2,  4,  -9 ],
          [ -6, -27, -8,  -6 ],
          [ -4,   8, 19,   6 ],
          [  4,  -5, -3, -13 ] ]
    b = [ -9, -76, -79, -70 ]
    print_system(a, b)
    x = seidel(a, b)
    for i in range(1, len(x)+1):
        print(f'x{get_script(str(i), "sub")} = %g' % round(x[i-1], 3))

print('1. Методом Гаусса решить СЛАУ')
print('2. Методом прогонки решить СЛАУ')
print('3. Методом простых итераций решить СЛАУ')
print('4. Методом Зейделя решить СЛАУ')
print('Пункт: ', end='')
[task1, task2, task3, task4][int(input())-1]()

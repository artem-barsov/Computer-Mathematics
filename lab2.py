from polynomial_algos import *

def task1():
    a = read_polynom()
    print('Производная:', polynom2str(derivat(a)))

def task2():
    a = read_polynom()
    roots = find_roots(a, chord, 1)
    print('Корни:', ', '.join(map(lambda x: '%g'%x, roots)))
    print('Наибольший корень:', '%g'%max(roots))

print('1. Найти производную полинома')
print('2. Методом хорд найти наибольший из корней уравнения')
print('Пункт: ', end='')
[task1, task2][int(input())-1]()

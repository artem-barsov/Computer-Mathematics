from polynomial_algos import *

def task1():
    a = read_polynom()
    print('Число а:', end=' ')
    t = float(input())
    print('Ответ:', f_x(a, t))

def task2():
    a = read_polynom()
    print('Число а:', end=' ')
    t = float(input())
    print('Частное:', polynom2str(quotient(a, t)))

def task3():
    a = read_polynom()
    print('Число а:', end=' ')
    t = float(input())
    print('Замена:', polynom2str(changed_vars(a, t), 'y'))

def task4():
    a = read_polynom()
    lower_bound, upper_bound = root_bounds(a)
    print('Нижняя граница корней: %g' % lower_bound)
    print('Верхняя граница корней: %g' % upper_bound)

def task5():
    a = read_polynom()
    print('Корни:', ', '.join(map(lambda x: '%g'%x, find_roots(a, dichotomy))))

print('1. Значение полинома в точке по схеме Горнера')
print('2. Частное от деления полинома на (x-a)')
print('3. Замена переменных в полиноме (x = y + a)')
print('4. Границы корней полинома')
print('5. Действительные корни полинома методом дихотомии')
print('Пункт: ', end='')
[task1, task2, task3, task4, task5][int(input())-1]()

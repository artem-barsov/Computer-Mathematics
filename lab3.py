from polynomial_algos import *

def task1():
    a = read_polynom()
    r, R = ring_bounds(a)
    print('Нижняя граница: %g' % r)
    print('Верхняя граница: %g' % R)

def task2():
    a = read_polynom()
    print('Верхняя граница: %g' % Lagrange_upper(a))

def task3():
    a = read_polynom()
    print('Верхняя граница: %g' % Newton_upper(a))

def task4():
    print('Первый полином:')
    a = read_polynom()
    print('Второй полином:')
    b = read_polynom()
    quot, rem = div_polys(a, b)
    print(f'({polynom2str(a)}) / ({polynom2str(b)}) = ...')
    print('Частное: ', polynom2str(quot))
    print('Остаток: ', polynom2str(rem))

print('1. Определить нижнюю и верхнюю границу кольца, внутри которого расположены все корни полинома')
print('2. Определить верхнюю границу положительных действительных корней полинома по методу Лагранжа')
print('3. Определение верхней границы положительных действительных корней полинома по методу Ньютона')
print('4. Реализовать деление полинома на полином произвольной степени (синтетическое деление)')
print('Пункт: ', end='')
[task1, task2, task3, task4][int(input())-1]()

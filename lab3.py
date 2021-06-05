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

def task5():
    a = read_polynom()
    print('Количество действительных корней:', cnt_real_roots(a))

def task6():
    e = 8
    a = read_polynom()
    l, R = root_bounds(a)
    print('Количество действительных корней:', cnt_real_roots(a))
    r = l+1
    while l <= R:
        cnt = cnt_real_roots(a, l, r)
        if cnt <= 1:
            if cnt == 1: print(f'({round(l, e)} ; {round(r, e)}):\t{dichotomy(a, l, r)}')
            l, r = r, r+1
        else: r = (l + r) / 2

print('1. Определить нижнюю и верхнюю границу кольца, внутри которого расположены все корни полинома')
print('2. Определить верхнюю границу положительных действительных корней полинома по методу Лагранжа')
print('3. Определение верхней границы положительных действительных корней полинома по методу Ньютона')
print('4. Реализовать деление полинома на полином произвольной степени (синтетическое деление)')
print('5. Определить количество действительных корней по методу Штурма')
print('6. Отделить корни полинома с помощью метода Штурма и вычислить их')
print('Пункт: ', end='')
[task1, task2, task3, task4, task5, task6][int(input())-1]()

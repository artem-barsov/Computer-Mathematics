from polynomial_algos import *
from matrix_algos import *
import system_algos

def task1():
    a = read_polynom()
    print('Производная:', polynom2str(d_dx(a)))

def task2():
    a = read_polynom()
    roots = find_roots(a, chord, 1)
    if roots:
        print('Корни:', ', '.join(map(lambda x: '%g'%x, roots)))
        print('Наибольший корень:', '%g'%max(roots))
    else:
        print('Невозможно найти корни методом хорд')

def task3():
    a = read_polynom()
    roots = find_roots(a, newton, 0.5)
    if roots:
        print('Корни:', ', '.join(map(lambda x: '%g'%x, roots)))
        print('Наибольший корень:', '%g'%max(roots))
    else:
        print('Невозможно найти корни методом Ньютона')

def task4():
    print('Матрица коэффициентов:')
    a = read_matrix()
    print('Полином:', matrix2polynomstr(a))
    print('Точка (x0,y0): ', end='')
    x, y = map(float, input().strip().split(' '))
    print('Ответ:', '%g'%f_x_y(a, x, y))

def task5():
    print('Матрица коэффициентов:')
    a = read_matrix()
    print('Полином:', matrix2polynomstr(a))
    print('Матрица производной по x:')
    da_dx = dmatrix_dx(a)
    print_matrix(da_dx)
    print('Производная по x:', matrix2polynomstr(da_dx))
    print('Матрица производной по y:')
    da_dy = dmatrix_dy(a)
    print_matrix(da_dy)
    print('Производная по y:', matrix2polynomstr(da_dy))

def task6():
    print('Действительная матрица коэффициентов:')
    real = read_matrix()
    print('Полином:', matrix2polynomstr(a))
    print('Мнимая матрица коэффициентов:')
    imag = read_matrix()
    print('Полином:', matrix2polynomstr(a))
    print('Ответ:', system_algos.newton(real, imag))

print('1. Найти производную полинома')
print('2. Методом хорд найти наибольший из корней уравнения')
print('3. Методом Ньютона найти наибольший из корней уравнения')
print('4. Найти значение полинома от двух переменных в точке по схеме Горнера')
print('5. Найти частные производные полинома от двух переменных')
print('6. Найти комплексный корень уравнения методом Ньютона для систем уравнений с двумя неизвестными')
print('Пункт: ', end='')
[task1, task2, task3, task4, task5, task6][int(input())-1]()

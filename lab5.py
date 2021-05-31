from polynomial_algos import *
from matrix_algos import *

def task1():
    print('Матрица:')
    a = read_matrix()
    print(polynom2str(char_polyn_LeVerrier(a), 'λ'))

def task2():
    print('Матрица:')
    a = read_matrix()
    max_eigenval, max_eigenvec = eigenval8vec(a)
    print('Максимальное по модулю собственное значение λ =', max_eigenval)
    print('Собственный вектор максимального по модулю собственного значения:')
    print_matrix(max_eigenvec)
    rest_eig = rest_eigenvalues3x3(a, max_eigenval)
    if not isinstance(rest_eig, tuple):
        print(f'Оставшееся собственное значение λ = {rest_eig} кратностью 2')
    elif not isinstance(rest_eig[0], tuple):
        print('Остальные собственные значения: {0}, {1}'.format(*rest_eig))
    else:
        z1, z2 = rest_eig
        print('Остальные собственные значения: ({0}+{1}i), ({2}+{3}i)'.format(*z1, *z2))

print('1. Нахождение характеристического полинома произвольной квадратной матрицы по методу Леверрье')
print('2. Нахождение собственных значений матрицы')
print('Пункт: ', end='')
[task1, task2][int(input())-1]()

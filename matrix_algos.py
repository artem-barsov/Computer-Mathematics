from polynomial_algos import d_dx

def read_matrix():
    print('Матрица коэффициентов:')
    a = []
    while (s:=input().strip())[-1] != ';':
        a.append([float(x) for x in s.split(' ')])
    a.append([float(x) for x in s[:-1].split(' ')])
    print('Полином:', matrix2polynomstr(a))
    return a

def matrix2polynomstr(a, x = 'x', y = 'y'):
    s = ''
    first = True
    for i in range(len(a)):
        degY = len(a) - i - 1
        for j in range(len(a[0])):
            if a[i][j] == 0: continue
            degX = len(a[0]) - j - 1
            if not first:
                s += ('+ ' if a[i][j] > 0 else '- ')
            elif a[i][j] < 0:
                s += '-'
            if abs(a[i][j]) != 1 or (degX==0 and degY==0):
                s += '%g'%abs(a[i][j])
            if degX > 0:
                s += x
                if degX != 1:
                    s += '^' + str(degX)
                if degY == 0:
                    s += ' '
            if degY > 0:
                if degX > 0: s += '*'
                s += y
                if degY != 1:
                    s += '^' + str(degY)
                s += ' '
            first = False
    return s

def dmatrix_dx(a):
    return list(map(d_dx, a))

def dmatrix_dy(a):
    return list(zip(*map(d_dx, zip(*a))))

def print_matrix(a):
    for ai in a:
        for aij in ai:
            print('%g'%aij, end=' ')
        print()

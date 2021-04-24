from polynomial_algos import d_dx, f_x_y

def read_matrix():
    print('Матрица коэффициентов:')
    a = []
    while (s:=input().strip())[-1] != ';':
        a.append([float(x) for x in s.split(' ')])
    a.append([float(x) for x in s[:-1].split(' ')])
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

def det2(a):
    return a[0][0]*a[1][1] - a[0][1]*a[1][0]

def Jacobian_matrix2(f, g, x, y):
    return [ [ f_x_y(dmatrix_dx(f), x, y), f_x_y(dmatrix_dy(f), x, y) ]
            ,[ f_x_y(dmatrix_dx(g), x, y), f_x_y(dmatrix_dy(g), x, y) ] ]

def mat_sum(a, b):
    if not isinstance(a, list): return a+b
    return list(map(mat_sum, a, b))

def mat_mul(a, b):
    if any([len(a)==1, len(b)==1, len(b[0])==1]):
        return [[sum([a[i][j]*b[j][k]
                    for j in range(len(a[0]))])
                for k in range(len(b[0]))]
            for i in range(len(a))]
    h1, c, v2 = len(a)//2, len(b)//2, len(b[0])//2
    A = [ [ [u[:c] for u in a[:h1]], [u[c:] for u in a[:h1]] ]
        , [ [u[:c] for u in a[h1:]], [u[c:] for u in a[h1:]] ] ]
    B = [ [ [u[:v2] for u in b[:c]], [u[v2:] for u in b[:c]] ]
        , [ [u[:v2] for u in b[c:]], [u[v2:] for u in b[c:]] ] ]
    return list(map(lambda a,b:a+b,
        mat_sum(mat_mul(A[0][0], B[0][0]), mat_mul(A[0][1], B[1][0])),
        mat_sum(mat_mul(A[0][0], B[0][1]), mat_mul(A[0][1], B[1][1]))
        )) + list(map(lambda a,b:a+b,
        mat_sum(mat_mul(A[1][0], B[0][0]), mat_mul(A[1][1], B[1][0])),
        mat_sum(mat_mul(A[1][0], B[0][1]), mat_mul(A[1][1], B[1][1]))
        ))

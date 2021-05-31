from polynomial_algos import d_dx, f_x_y, quotient
from math import log10, sqrt

def read_matrix():
    a = []
    while (s:=input().strip())[-1] != ';':
        a.append([float(x) for x in s.split(' ')])
    a.append([float(x) for x in s[:-1].split(' ')])
    return a

def matrix2polynomstr(a, x = 'x', y = 'y'):
    def get_super(s):
        return s.translate(s.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹'))
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
                    s += get_super(str(degX))
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

def print_matrix(a, e = 3):
    for ai in a:
        for aij in ai:
            print('%g'%(round(aij, e)+0.), end=' ')
        print()

def Jacobian_matrix2(f, g, x, y):
    return [ [ f_x_y(dmatrix_dx(f), x, y), f_x_y(dmatrix_dy(f), x, y) ]
            ,[ f_x_y(dmatrix_dx(g), x, y), f_x_y(dmatrix_dy(g), x, y) ] ]

def quarters2matrix(a):
    return list(map(lambda a,b:a+b, a[0][0], a[0][1])) + list(map(lambda a,b:a+b, a[1][0], a[1][1]))

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
    return quarters2matrix( [ [ mat_sum(mat_mul(A[0][0], B[0][0]), mat_mul(A[0][1], B[1][0])),
                                mat_sum(mat_mul(A[0][0], B[0][1]), mat_mul(A[0][1], B[1][1])) ],
                              [ mat_sum(mat_mul(A[1][0], B[0][0]), mat_mul(A[1][1], B[1][0])),
                                mat_sum(mat_mul(A[1][0], B[0][1]), mat_mul(A[1][1], B[1][1])) ] ] )

def inverse_matrix(a):
    if len(a) != len(a[0]):
        raise Exception('Невозможно найти обратную матрицу для неквадратной матрицы')
    if len(a) == 1:
        if a[0][0] == 0:
            raise Exception('Невозможно найти обратную матрицу: один из промежуточных блоков особенный')
        return [[1/a[0][0]]]
    S = [[a[0][0]]]                                       # k x k
    S_inv = inverse_matrix(S)                             # k x k
    for k in range(len(a)-1):
        a12 = [[ai[len(S)]] for ai in a[:len(S)]]         # k x 1
        a21 = [a[len(S)][:len(S)]]                        # 1 x k
        a22 = [[a[len(S)][len(S)]]]                       # 1 x 1
        X = mat_mul(S_inv, a12)                           # k x 1
        Y = mat_mul(a21, S_inv)                           # 1 x k
        if a22[0][0] == mat_mul(a21, X)[0][0]:
            raise Exception('Невозможно найти обратную матрицу: один из промежуточных блоков особенный')
        T_inv = [[1/(a22[0][0] - mat_mul(a21, X)[0][0])]] # 1 x 1
        XT_inv = mat_mul(X, T_inv)                        # k x 1
        T_invY = mat_mul(T_inv, Y)                        # 1 x k
        XT_invY = mat_mul(XT_inv, Y)                      # k x k
        for i in range(len(XT_inv)):
            for j in range(len(XT_inv[0])):
                XT_inv[i][j] = -XT_inv[i][j]
        for i in range(len(T_invY)):
            for j in range(len(T_invY[0])):
                T_invY[i][j] = -T_invY[i][j]
        S_inv = quarters2matrix([ [ mat_sum(S_inv, XT_invY), XT_inv ]     # k x k, k x 1
                                , [ T_invY,                  T_inv  ] ] ) # 1 x k, 1 x 1
        S = quarters2matrix([ [ S,   a12 ]
                            , [ a21, a22 ] ] )
    return S_inv

def nonzero_el(a, i0 = 0, j0 = 0):
    for i in range(i0, len(a)):
        for j in range(j0, len(a[0])):
            if a[i][j] != 0:
                return i, j
    return None

def det(a):
    if len(a) != len(a[0]):
        raise Exception('Невозможно найти определитель неквадратной матрицы')
    ret = 1
    for k in range(len(a)):
        if not (idx := nonzero_el(a, k, k)):
            return 0
        i0, j0 = idx
        if i0 != k:
            a[k], a[i0] = a[i0], a[k]
            ret = -ret
        if j0 != k:
            for ai in a:
                ai[k], ai[j0] = ai[j0], ai[k]
            ret = -ret
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a[0])):
                a[i][j] -= a[i][k]*a[k][j]/a[k][k]
        ret *= a[k][k]
    return ret

def track(a):
    return sum([a[i][i] for i in range(len(a))])

def char_polyn_LeVerrier(a):
    ak = a
    s = [track(a)]
    for _ in range(len(a)-1):
        ak = mat_mul(ak,a)
        s.append(track(ak))
    p = []
    for i in range(len(a)):
        p.append(-(s[i] + sum([p[j]*s[i-j-1] for j in range(i)])) / (i+1))
    return [1] + p

def normalized_v(v):
    norm = sqrt(sum([x**2 for x in v]))
    return [vx/norm for vx in v]

def eigenval8vec(a, e = 0.001):
    y = [[1] for _ in range(len(a))]
    ak = a
    yk_prev = y
    eigenval_prev = 0
    eigenval = e
    while abs(eigenval - eigenval_prev) >= e:
        eigenval_prev = eigenval
        yk = mat_mul(ak, y)
        eigenval = sum([yk[i][0]/yk_prev[i][0] for i in range(len(a))]) / len(a)
        yk_prev = yk
        ak = mat_mul(ak, a)
    eigenvec = list(zip(normalized_v(list(*zip(*yk)))))
    return round(eigenval, int(log10(1/e))), eigenvec

def rest_eigenvalues3x3(m, max_eigenval, e = 0.001):
    a,b,c = quotient(char_polyn_LeVerrier(m), max_eigenval)
    D = b**2 - 4*a*c
    if D == 0:
        return round(-b/(2*a), int(log10(1/e)))
    if D > 0:
        return round((-b + sqrt(D))/(2*a), int(log10(1/e))), round((-b - sqrt(D))/(2*a), int(log10(1/e)))
    if D < 0:
        return (round(-b/(2*a), int(log10(1/e))), round(sqrt(-D)/(2*a)), int(log10(1/e))), \
        (round(-b/(2*a), int(log10(1/e))), round(-sqrt(-D)/(2*a), int(log10(1/e))))

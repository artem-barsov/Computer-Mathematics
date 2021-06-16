from matrix_algos import Jacobian_matrix2, det, dmatrix_dx, dmatrix_dy, mat_mul, mat_sum
from polynomial_algos import f_x_y, polynom2str
from math import log10

def newton(f, g, e = 0.001):
    x0, y0 = 1.2, 1.7 # подбор
    J = det(Jacobian_matrix2(f, g, x0, y0))
    D = [ [ f_x_y(f, x0, y0), f_x_y(dmatrix_dy(f), x0, y0) ]
        , [ f_x_y(g, x0, y0), f_x_y(dmatrix_dy(g), x0, y0) ] ]
    xn = x0 - det(D) / J
    D = [ [ f_x_y(dmatrix_dx(f), x0, y0), f_x_y(f, x0, y0) ]
        , [ f_x_y(dmatrix_dx(g), x0, y0), f_x_y(g, x0, y0) ] ]
    yn = y0 - det(D) / J
    while abs(xn - x0) + abs(yn - y0) >= e:
        J = det(Jacobian_matrix2(f, g, xn, yn))
        D1= [ [ f_x_y(f, xn, yn), f_x_y(dmatrix_dy(f), xn, yn) ]
            , [ f_x_y(g, xn, yn), f_x_y(dmatrix_dy(g), xn, yn) ] ]
        D2= [ [ f_x_y(dmatrix_dx(f), xn, yn), f_x_y(f, xn, yn) ]
            , [ f_x_y(dmatrix_dx(g), xn, yn), f_x_y(g, xn, yn) ] ]
        x0, xn = xn, xn-det(D1) / J
        y0, yn = yn, yn-det(D2) / J
    return round(xn, int(log10(1/e))), round(yn, int(log10(1/e)))

def print_system(a, b):
    for i in range(len(a)):
        print('{', polynom2str(a[i], script='sub'), '=', b[i])

def Gauss(a, b):
    for i in range(len(a)-1):
        r = a[i][i]
        for j in range(len(a)-i-1):
            t = a[i+j+1][i]
            for k in range(i, len(a)):
                a[i+j+1][k] -= a[i][k]/r*t
            b[i+j+1] -= b[i]/r*t
    x = []
    for i in range(len(a)-1, -1, -1):
        x.append((b[i]-sum([a[i][-j-1]*x[j] for j in range(len(x))]))/a[i][-len(x)-1])
    return x[::-1]

def tridiagonal_algo(a, b):
    alfa = [0] + [a[i+1][i] for i in range(len(a)-1)]
    beta = [-a[i][i] for i in range(len(a))]
    gamma = [a[i][i+1] for i in range(len(a)-1)] + [0]
    P = [gamma[0] / beta[0]]
    Q = [-b[0] / beta[0]]
    for i in range(1, len(a)-1):
        P.append(gamma[i] / (beta[i] - alfa[i] * P[i-1]))
        Q.append((alfa[i] * Q[i-1] - b[i]) / (beta[i] - alfa[i] * P[i-1]))
    x = [(alfa[-1] * Q[-1] - b[-1]) / (beta[-1] - alfa[-1] * P[-1])]
    for i in range(len(a)-2, -1, -1):
        x.append(P[i] * x[-1] + Q[i])
    return x[::-1]

def simple_iteration(a, b, e = 0.01):
    b = [[b[i] / a[i][i]] for i in range(len(b))]
    a = [ [-a[i][j]/a[i][i] for j in range(len(a))] for i in range(len(a)) ]
    for i in range(len(a)): a[i][i] = 0
    x = b.copy()
    while True:
        x_prev = x.copy()
        x = mat_sum(mat_mul(a, x), b)
        if (sum([abs(x[i][0]-x_prev[i][0]) for i in range(len(x))]) < e):
            break
    return list(*zip(*x))

def seidel(a, b, e = 0.01):
    b = [[b[i] / a[i][i]] for i in range(len(b))]
    a = [ [-a[i][j]/a[i][i] for j in range(len(a))] for i in range(len(a)) ]
    for i in range(len(a)): a[i][i] = 0
    x = b.copy()
    while True:
        x_prev = x.copy()
        x = []
        for i in range(len(a)):
            xi = 0
            # print('len(x) =', len(x))
            for j in range(len(x)):
                # print(f'xi += a[{i}][{j}] * x[{j}][{0}]')
                xi += a[i][j] * x[j][0]
            for j in range(len(x), len(a)):
                xi += a[i][j] * x_prev[j][0]
            x.append([xi + b[i][0]])
        # x = mat_sum(mat_mul(a, x), b)
        if (sum([abs(x[i][0]-x_prev[i][0]) for i in range(len(x))]) < e):
            break
    return list(*zip(*x))

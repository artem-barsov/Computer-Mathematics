from math import ceil, log10

def f(a, t):
    b = 0
    for ai in a:
        b = ai + b*t
    return b

def f_x_y(a, x, y):
    return f(map(lambda a: f(a, x), a), y)

def polynom2str(a, x = 'x'):
    def get_super(s):
        return s.translate(s.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹'))
    s = ''
    for i in range(len(a)):
        if a[i] == 0: continue
        deg = len(a) - i - 1
        if i != 0:
            s += ('+ ' if a[i] > 0 else '- ')
        elif a[i] < 0:
            s += '-'
        if abs(a[i]) != 1 or deg == 0:
            s += '%g'%abs(a[i])
        if deg > 0:
            s += x
            if deg != 1:
                s += get_super(str(deg))
            s += ' '
    return s if s!='' else '0'

def read_polynom():
    print('Коэффициенты: ')
    a = [float(x) for x in input().strip().split(' ')]
    while a[0]==0: del a[0]
    print('Степень полинома =', len(a)-1)
    print('Полином:', polynom2str(a))
    return a

def quotient(a, t):
    ret = [a[0]]
    for ai in a[1:-1]:
        ret.append(ai + ret[-1]*t)
    return ret

def changed_vars(a, t):
    ret = []
    for i in range(len(a)):
        for j in range(1, len(a)-i):
            a[j] += a[j-1]*t
        ret.append(a[-i-1])
    return ret[::-1]

def div_polys(a, b):
    a = [ai/b[0] for ai in a]
    b = [-bi/b[0] for bi in b[1:]]
    ret = [0 for _ in range(len(b))]
    for i in range(len(a)):
        ret.append(
            a[i]+sum([b[j]*ret[-(j+1)]
                for j in range( len(b)-1,
                    -1+(i-len(a)+len(b))*(i>=len(a)-len(b)+1),
                    -1)]))
    return ret[len(b):-len(b)], ret[-len(b):]

def Lagrange_upper(a):
    if a[0]<0: a = [-ai for ai in a]
    neg_ids = [i for i,x in enumerate(a) if x<0]
    A = -min([a[i] for i in neg_ids])
    return 1 + pow(A/a[0], 1/neg_ids[0])

def Newton_upper(a, step = 1):
    if a[0]<0: a = [-ai for ai in a]
    dn_dx = [a]
    for _ in range(len(a)-3):
        dn_dx.append(d_dx(dn_dx[-1]))
    c = 0
    for dk_dx in dn_dx[::-1]:
        while f(dk_dx, c) < 0: c += step
    return c

def root_bounds(a):
    def find_upper(a, step = 1):
        upp = ceil(Lagrange_upper(a))
        while all([x>=0 for x in quotient(a, upp-step)]):
            upp -= step
        return upp
    a = [x/a[0] for x in a]
    b = [a[i]*(1 if (len(a)-i-1)%2==0 else -1) for i in range(len(a))]
    b = [x/b[0] for x in b]
    if any([x<0 for x in a]): # Has positive roots
        if any([x<0 for x in b]): # Has negative roots
            return -find_upper(b), find_upper(a)
        else: # No negative roots
            return (1 / find_upper([x/a[-1] for x in a[::-1]])), find_upper(a)
    else: # No positive roots
        return -find_upper(b), -1 / find_upper([x/b[-1] for x in b[::-1]])

def sign(x): return (x > 0) - (x < 0)

def d_dx(a):
    return [ai*(len(a)-i-1) for (i, ai) in enumerate(a[:-1])]

def dichotomy(a, l, r, e = 0.001):
    if sign(f(a, l)) == sign(f(a, r)):
        return None
    while r - l >= e:
        m = (l + r) / 2
        if sign(f(a, l)) != sign(f(a, m)):
            r = m
        else:
            l = m
    return round((l + r) / 2, int(log10(1/e)))

def chord(a, l, r, e = 0.001):
    if sign(f(a, l)) == sign(f(a, r)):
        return None
    while r - l >= e:
        m = l - f(a, l) * (r - l) / (f(a, r) - f(a, l))
        f_m = f(a, m)
        if abs(f_m) < e: return round(m, int(log10(1/e)))
        if sign(f(a, l)) != sign(f_m):
            r = m
        else:
            l = m
    return round(l, int(log10(1/e)))

def newton(a, l, r, e = 0.001):
    da_dx = d_dx(a)
    d2a_dx2 = d_dx(da_dx)
    if sign(f(a, l)) == sign(f(a, r)):
        return None
    if sign(f(da_dx, l)) != sign(f(da_dx, r)):
        return None
    if sign(f(d2a_dx2, l)) != sign(f(d2a_dx2, r)):
        return None
    x0 = l if sign(f(a, l)) == sign(f(da_dx, l)) else r
    if f(da_dx, x0) == 0: return round(x0, int(log10(1/e)))
    xn = x0 - f(a, x0)/f(da_dx, x0)
    while abs(xn - x0) >= e:
        x0, xn = xn, xn-f(a, xn)/f(da_dx, xn)
    return round(xn, int(log10(1/e)))

def find_roots(a, method = dichotomy, step = 0.3):
    l, r = root_bounds(a)
    ret = []
    while l <= r:
        x = method(a, l, l+step)
        if x is not None: ret.append(x)
        l += step
    return sorted(set(ret))

def ring_bounds(a):
    A = max([abs(ai) for ai in a[1:]])
    R = 1 + A / abs(a[0])
    B = max([abs(ai) for ai in a[:-1]])
    i = 1
    while a[-i] == 0: i += 1
    r = 1 / (1 + B / abs(a[-i]))
    return r, R

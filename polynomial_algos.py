from math import ceil, log10

def f_x(a, t):
    b = 0
    for ai in a:
        b = ai + b*t
    return b

def polynom2str(a, x = 'x'):
    s = ''
    for i in range(len(a)):
        deg = len(a)-i-1
        if a[i] == 0: continue
        if i != 0:
            s += ('+ ' if a[i] > 0 else '- ')
        elif a[i]<0:
            s += '-'
        if abs(a[i]) != 1 or deg == 0:
            s += '%g'%abs(a[i])
        if deg > 0:
            s += x
            if deg != 1:
                s += '^' + str(deg)
            s += ' '
    return s

def read_polynom():
    print('Коэффициенты: ')
    a = [float(x) for x in input().strip().split(' ')]
    while a[0]==0: del a[0]
    print('Степень полинома =', len(a)-1)
    print('Полином:', polynom2str(a))
    return a

def quotient(a, t):
    ret = [a[0]]
    for ai in a[1:]:
        ret.append(ai + ret[-1]*t)
    return ret

def changed_vars(a, t):
    ret = []
    for i in range(len(a)):
        for j in range(1, len(a)-i):
            a[j] += a[j-1]*t
        ret.append(a[-i-1])
    return ret[::-1]

def root_bounds(a):
    def find_upper(a, step = 1):
        neg_ids = [i for i,x in enumerate(a) if x<0]
        A = -min([a[i] for i in neg_ids])
        upp = ceil(1 + pow(A/a[0], 1/neg_ids[0]))
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

def dichotomy(a, l, r, e = 0.001):
    if sign(f_x(a, l)) == sign(f_x(a, r)):
        return None
    while r - l >= e:
        m = (l + r) / 2
        if sign(f_x(a, l)) != sign(f_x(a, m)):
            r = m
        else:
            l = m
    return round((l + r) / 2, int(log10(1/e)))

def chord(a, l, r, e = 0.001):
    if sign(f_x(a, l)) == sign(f_x(a, r)):
        return None
    while r - l >= e:
        m = l - f_x(a, l) * (r - l) / (f_x(a, r) - f_x(a, l))
        if abs(f_m:=f_x(a, m)) < e: return round(m, int(log10(1/e)))
        if abs(f_l:=f_x(a, l)) < e: return round(l, int(log10(1/e)))
        if abs(f_r:=f_x(a, r)) < e: return round(r, int(log10(1/e)))
        if sign(f_x(a, l)) != sign(f_m):
            r = m
        else:
            l = m
    return round(l, int(log10(1/e)))

def find_roots(a, method = dichotomy, step = 0.3):
    l, r = root_bounds(a)
    ret = []
    while l < r:
        x = method(a, l, l+step)
        if x is not None: ret.append(x)
        l += step
    return ret

def derivat(a):
    return [ai*(len(a)-i-1) for (i, ai) in enumerate(a[:-1])]

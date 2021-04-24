from matrix_algos import Jacobian_matrix2, det2, dmatrix_dx, dmatrix_dy
from polynomial_algos import f_x_y
from math import log10

def newton(f, g, e = 0.001):
    x0, y0 = 1.2, 1.7 # подбор
    J = det2(Jacobian_matrix2(f, g, x0, y0))
    D = [ [ f_x_y(f, x0, y0), f_x_y(dmatrix_dy(f), x0, y0) ]
        , [ f_x_y(g, x0, y0), f_x_y(dmatrix_dy(g), x0, y0) ] ]
    xn = x0 - det2(D) / J
    D = [ [ f_x_y(dmatrix_dx(f), x0, y0), f_x_y(f, x0, y0) ]
        , [ f_x_y(dmatrix_dx(g), x0, y0), f_x_y(g, x0, y0) ] ]
    yn = y0 - det2(D) / J
    while abs(xn - x0) + abs(yn - y0) >= e:
        J = det2(Jacobian_matrix2(f, g, xn, yn))
        D1= [ [ f_x_y(f, xn, yn), f_x_y(dmatrix_dy(f), xn, yn) ]
            , [ f_x_y(g, xn, yn), f_x_y(dmatrix_dy(g), xn, yn) ] ]
        D2= [ [ f_x_y(dmatrix_dx(f), xn, yn), f_x_y(f, xn, yn) ]
            , [ f_x_y(dmatrix_dx(g), xn, yn), f_x_y(g, xn, yn) ] ]
        x0, xn = xn, xn-det2(D1) / J
        y0, yn = yn, yn-det2(D2) / J
    return round(xn, int(log10(1/e))), round(yn, int(log10(1/e)))

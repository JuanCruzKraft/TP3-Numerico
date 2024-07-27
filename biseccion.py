import matplotlib.pylab as plt
import numpy as np
from math import sqrt

def f(x):
    return sqrt(x)-15

def biseccion(a,b,maxit,tol) :
    i = 0
    err = []
    m = (a + b / 2.0)  # intervalo
    q = f(m) # residuo
    error = abs((m - a) / 2) # error
    #iterar para encontrar la solucion
    while(i < maxit and (not bool(err) or abs(err[-1]) > tol)) :
        m = (a + b) / 2
        q = f(m)
        if np.sign(f(a)) == np.sign(q):
            a = m
        else:
            b = m
        error = abs((a - b) / 2)
        err.append(error)
        i+=1
    return i, m, err


i, m, e = biseccion(100,300,50,0.001)
print(f"{i} iteraciones, m = {m}, error: {e}")
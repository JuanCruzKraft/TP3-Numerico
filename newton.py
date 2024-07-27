import math

# x[i+1] = x[i] - f(x[i])/f'(x[i])

# derivada de f
def g(x):
    return -math.exp(-x) - 1
def f(x):
    return math.exp(-x) - x
 
def newton(x0, maxit, tol):
    i = 0
    error = []
    x = x0
    y = f(x0)
    error.append(y)
    while(i < maxit and abs(error[-1])>tol):
        x -= f(x)/g(x)
        y = f(x)
        error.append(y)
        i+=1
    return x, error

#Xn-1 = Xt
def newton_secante(x0, x1, maxit, tol, func, h):
    error = []
    Xn = x1
    Xt = x0
    Y = func(x0,h)
    error.append(Y)
    i = 0
    while(i<maxit):
        if(abs(error[-1])>tol):
            break
        aux = Xn
        Xn = Xn - ((Xn - Xt)/(func(Xn,h)-func(Xt,h)))*func(Xn,h)
        Y = func(Xn,h)
        Xt = aux
        error.append(Y)
        i+=1
        print(i, Xn, Xt)
    return Xn


def newton_secante2(f, df, d2f, x0, a_n, b_n, tol=1e-7, max_iter=100):
    x_n = x0
    history = [x_n]
    
    for i in range(max_iter):
        P_NR = df(x_n)
        P_sec = (f(b_n) - f(a_n)) / (b_n - a_n)
        
        w_a = d2f(a_n) / (d2f(a_n) + d2f(b_n))
        w_b = d2f(b_n) / (d2f(a_n) + d2f(b_n))
        
        MNS = w_a * P_NR + w_b * P_sec
        
        x_n1 = x_n - f(x_n) / MNS
        history.append(x_n1)
        
        if abs(f(x_n1)) < tol:
            break
            
        x_n = x_n1
    
    print(f"Cantidad de iteraciones: {len(history) - 1}")
    return x_n, history

    

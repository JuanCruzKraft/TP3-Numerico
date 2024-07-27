def regula_falsi(f, a, b, tol=1e-6, max_iter=100):
    if f(a) * f(b) >= 0:
        raise ValueError("Inicio y fin deben tener distinto signo.")
    
    for iteration in range(max_iter):
        # Calculate the point of intersection with the x-axis
        c = b - (f(b) * (b - a)) / (f(b) - f(a))
        
        # Check if the result is within the tolerance
        if abs(f(c)) < tol:
            return (c, iteration)
        
        # Update the interval
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    
    raise RuntimeError("Máximo número de iteraciones.")

def newton_raphson(f, f_derivada, x0=0, tol=1e-6, max_iter=1000):
    x = x0
    for iteration in range(max_iter):
        fx = f(x)
        dfx = f_derivada(x)
        
        if dfx == 0:
            raise ValueError("The derivative is zero at x = {}. No solution found.".format(x))
        
        x_new = x - fx / dfx
        
        if abs(x_new - x) < tol:
            return (x_new, iteration)
        
        x = x_new
    
    raise RuntimeError("Maximum number of iterations reached. No solution found.")

# Implementación del método Newton-Secante con segunda derivada
def newton_secante(f, f_derivada, f_segunda_derivada, x0, x1, tol=1e-6, max_iter=1000):
    for i in range(max_iter):
        f0, f1 = f(x0), f(x1)
        df1 = f_derivada(x1)
        ddf0, ddf1 = f_segunda_derivada(x0), f_segunda_derivada(x1)
        alpha = ddf1 / (ddf1 + ddf0)
        mns = (1 - alpha) * df1 + alpha * (f1 - f0) / (x1 - x0)
        x2 = x1 - f1 / mns
        if abs(x2 - x1) < tol:
            return x2, i+1
        x0, x1 = x1, x2
    raise ValueError("No convergió")

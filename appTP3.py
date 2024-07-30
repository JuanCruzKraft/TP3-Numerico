import metNumpy as mN
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from sympy.utilities.lambdify import lambdify


# Definimos una función de ejemplo y sus derivadas
def f(x):
    return x**3 - 6*x**2 + 11*x - 6  #funcion ejemplo

def f_derivada(x):
    return 3*x**2 - 12*x + 11

def f_segunda_derivada(x):
    return 6*x - 12


entrar = input('Ejecutar inciso a y b (y/n): ')
if entrar == 'y':
    # Parámetros iniciales
    x0, x1 = 0.0, 1.5  # Valores iniciales


    # Imprimimos raices 
    rootNS, iterationsNS, erroresNS = mN.newton_secante(f, f_derivada, f_segunda_derivada, x0, x1)
    print(f"Por Newton-Secante se encontró la raíz: {rootNS}, en {iterationsNS} iteraciones. ")
    print(f"El vector de errores que contiene cada error calculado en cada iteracion es el siguiente: {erroresNS}")

    rootRF, iterationRF, erroresRF = mN.regula_falsi(f, 0.0, 1.5)
    print(f"Por Regula-Falsi se encontró la raíz: {rootRF}, en {iterationRF} iteraciones. ")
    print(f"El vector de errores que contiene cada error calculado en cada iteracion es el siguiente: {erroresRF}")
    rootNR, iterationNR , erroresNR= mN.newton_raphson(f, f_derivada)
    print(f"Por Newton-Raphson se encontró la raíz: {rootNR}, en {iterationNR} iteraciones. ")
    print(f"El vector de errores que contiene cada error calculado en cada iteracion es el siguiente: { erroresNR}")

    # Visualización gráfica
    x = np.linspace(-10, 10, 400)  # Ajusta el rango de x para cubrir de -10 a 10
    y = f(x)

    # Configura límites para el eje Y
    y_min = -10
    y_max = 10

    plt.plot(x, y, label='f(x)')
    plt.axhline(0, color='black', lw=0.5)
    plt.axvline(rootNS, color='red', linestyle='--', label=f'Raíz aproximada: {rootNS:.4f}')

    # Graficar las iteraciones
    points_x = [x0, x1]
    points_y = [f(x0), f(x1)]

    for i in range(1, iterationsNS + 1):
        f0, f1 = f(points_x[-2]), f(points_x[-1])
        df1 = f_derivada(points_x[-1])
        ddf0, ddf1 = f_segunda_derivada(points_x[-2]), f_segunda_derivada(points_x[-1])
        alpha = ddf1 / (ddf1 + ddf0)
        mns = (1 - alpha) * df1 + alpha * (f1 - f0) / (points_x[-1] - points_x[-2])
        x_next = points_x[-1] - f1 / mns
        points_x.append(x_next)
        points_y.append(f(x_next))
        plt.plot(points_x[-3:], points_y[-3:], marker='o', linestyle='-', color='gray', alpha=0.5)

    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend()
    plt.title('Método Newton-Secante')
    plt.grid(True)

    plt.xticks(np.arange(-10, 11, 1))
    plt.yticks(np.arange(-10, 11, 1))

    # Ajusta los límites del eje Y para mantenerlos dentro del rango deseado
    plt.xlim(-10, 10)  # Ajuste del rango del eje X
    plt.ylim(y_min, y_max)  # Ajuste del rango del eje Y

    # plt.show()

####### Ejercicio C
print("---------------- INICIO EJERCICIO C ----------------")
def symbolic_derivative(func, symbol):

    symbol = sp.symbols(symbol)

    # Convert the Python function to a symbolic expression
    func_expr = func(symbol)
    
    # Calculate the derivative
    derivative_expr = sp.diff(func_expr, symbol)
    
    # Convert the symbolic derivative to a Python function
    derivative_func = lambdify(symbol, derivative_expr, 'numpy')
    
    return derivative_func


Ka_acetico = (10**(4.756))
Ka_amonico = 1/(10**(9.25))
Kw = 10**(-14)

entrada = input("desea graficar (y/n): ")
if entrada == 'y':
    # Rango de H para graficar
    H_values = np.logspace(-14, 0, 400)

    # Graficar f(H) para diferentes valores de CH3COOH y NH3
    for i in range(0, 10):
        CH3COOH = i + 1
        NH3 = 10 - i
        f_values = f(H_values)
        plt.plot(H_values, f_values, label=f'CH3COOH={CH3COOH}, NH3={NH3}')

    # Configuración de la gráfica
    plt.xscale('log')
    plt.xlabel('H')
    plt.ylabel('f(H)')
    plt.title('Gráfica de f(H) para diferentes valores de CH3COOH y NH3')
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.legend()
    plt.grid(True, which='both', ls='--', lw=0.5)
    plt.show()

resultados = []

for i in range(0,10):
    CH3COOH = i + 1
    NH4     = 10 - i

    CH3COOH /= 1000
    NH4 /= 1000

    def f(H):
        return H - (Kw/H) + ((NH4*H)/(Ka_amonico + H)) - (CH3COOH/((Ka_acetico * H) + 1))

    

    ## alternativa para derivar
    """ def df(H):
        h = 1e-6
        return (f(H + h) - f(H)) / h
    def d2f(H):
        h = 1e-6
        return (df(H + h) - df(H)) / h """

    ## Derivamos f
    f_derivada = symbolic_derivative(f, 'H')
    f_segunda_derivada = symbolic_derivative(f_derivada, 'H')

    """ print(f"KAacetico-KAamonico: {Ka_acetico} - {Ka_amonico} <-> Kw = {Kw}")
    print(f"CH3COOH = {CH3COOH} -- NH4 = {NH4}") """

    rootH, iterationsH , erroresH= mN.newton_secante(f, f_derivada, f_segunda_derivada, 1e-10, 2e-10, 1e-10) 
    pH = -np.log10(rootH)
    print(f"rootH: {rootH}; pH: {pH}; CH3COOH: {CH3COOH}, NH3: {NH4}, ITERACIONES: {iterationsH}")
    resultados.append(pH)
    
## Graficar resultados
i = [1,2,3,4,5,6,7,8,9,10]
plt.plot(i,resultados,'o-', label='pH en cada columna i')
plt.xlabel('i')
plt.ylabel('pH')
plt.title('Grafico de pH x columna')
plt.grid(True)
plt.show()
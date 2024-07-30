import numpy as np

## Funcion que determinara los valores de la matriz jacobiana
def resJ(df, x0, N, sum):
    # lleno la matriz Jacobiana de NxN con ceros
    J0 = np.zeros((N,N))
    
    # calcular para cada posicion el valor de la derivada correspondiente
    for i in range(N): #filas
        for j in range(-sum,sum+1): #columnas
            if(i+j>=0 and i+j <= N-1):
                J0[i,i+j] = df(x0[i+j], i+j, j)
    return J0

## Funcion que determinara los valores del vector de funciones F de tamaño N
def resF(f, x0, N):
    # llenar de ceros el vector F de tamaño N
    F0 = np.zeros(N)
    
    # calcular para cada posicion (i) el valor de f(H)
    for i in range(N):
        F0[i] = f(x0,i)
        
    return F0

#### metodo para resolver con Newton Jacobiano
## F (contiene las funciones para evaluar)
## J (contiene las funciones a evaluar)
## x0 (contiene el vector de condiciones iniciales)
def newtonJacobiano(f, df, x0, N, sum, iteraciones):
    i = 0
    for i in range(iteraciones):
        # evaluar J con x0 -> J(x0)
        J0 = resJ(df, x0 ,N, sum)
        
        # evaluar F con x0 -> F(x0)
        F0 = resF(f,x0,N)
        
        #resolver sistema Ax=b ->J=A b=F
        y0 = np.linalg.solve(J0, -F0).reshape(-1) # METODO DIRECTO de numpy para resolver Ax=b factorizacion LU
        
        #guardar vector de resolucion en x0 y volver a empezar.
        x0 += y0
    # retornar el resultado
    return x0, i

#### Metodo 2 de resolucion de sistema
#### Metodo 3 de resolucion de sistema

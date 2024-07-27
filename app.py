import newton as nw

# TP 3 - Metodo Newton Secante
KW = 10**-8
KA_ACETICO = 4.756
KA_AMONIO = 9.25

def conAmo(y):
    return 10 - y
def conAce(y):
    return 1 + y

def f(x, y):
    t0 = KW/x
    t1 = (conAmo(y)*x)/(10.0**KA_AMONIO)
    t2 = conAce(y)/((10.0**KA_ACETICO)*x)
    
    return x - t0 + t1 - t2

# --------------------------------------------------
res = []
for i in range(10):
    res.append(nw.newton_secante(1e-6,1e-7,100,1e-20, f, i))
print(res)
    


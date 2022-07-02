import numpy as np
import matplotlib.pyplot as plt  # Llamad de módulos necesarios


def sustA(A, b):  # Función de sustitución hacia adelante
    n = len(b)
    y = np.zeros([n, 1])
    for i in range(n):
        suma = 0
        for j in range(i):
            suma += A[i][j]*y[j][0]
        y[i][0] = (b[i][0]-suma)/A[i][i]
    return y  # Retorna el vector solución de la sustitución


def sustB(A, b):  # Función de sustitución hacia atrás
    n = len(b)
    x = np.zeros([n, 1])
    for i in reversed(range(n)):
        suma = 0
        for j in range(i+1, n):
            suma += A[i][j]*x[j][0]
        x[i][0] = (b[i][0]-suma)/A[i][i]

    return x  # Retorna el vector solución de la sustitución


def gaussSimple(A, b):  # Función que soluciona un sistema de ecuaciones por Gauss simple
    n = len(A)  # Longitud de la matriz
    for k in range(n):
        for i in range(k):
            factor = (A[k][i])/(A[i][i])
            for j in range(n):
                A[k][j] -= (factor*A[i][j])
            # Estos términos es la implementación de las fórmulas disponibles para hacer este método (pivoteo)
            b[k][0] = b[k][0] - (factor*b[i][0])

    x = sustB(A, b)  # Se llama a la función para hacer la sustitución
    return x


def crout(A, b):  # Función que soluciona un sistema de ecuaciones por descomposición L, U (Crout)
    n = len(A)  # Longitud de la matriz
    U = np.zeros([n, n])
    # Se crean dos matrices de zeros con el tamaño necesario
    L = np.zeros([n, n])

    for j in range(n):
        U[j][j] = 1
        for i in range(j, n):
            sum = 0
            for k in range(j):
                sum += L[i][k]*U[k][j]
            L[i][j] = A[i][j]-sum
        for k in range(j+1, n):
            sum = 0
            for i in range(j):
                sum += L[j][i]*U[i][k]
            U[j][k] = (A[j][k]-sum)/L[j][j]

    # Como al utilizar este método se necesitan hacer dos sustituciones, se hacen en esta unica linea con las funciones ya hechas
    x = sustB(U, sustA(L, b))

    return x


# Para este método de interpolación se piden los datos iniciales y se dan por defectos las iteraciones y tolerancia
def interpolacion(f, x0, x1, x2, n=100, tol=10**-8):
    for i in range(n):
        v1 = f(x0)*(x1**2-x2**2)+f(x1)*(x2**2-x0**2)+f(x2)*(x0**2-x1**2)
        v2 = 2*f(x0)*(x1-x2)+2*f(x1)*(x2-x0)+2*f(x2)*(x0-x1)
        x3 = v1/v2
        if x3 > x1:
            if f(x3) > f(x1):
                x0 = x1
                x1 = x3
            else:
                x2 = x3
        else:
            if f(x3) > f(x1):
                x2 = x1
                x1 = x3

            if f(x3) < f(x1):
                x0 = x3
        v1p = f(x0)*(x1**2-x2**2)+f(x1)*(x2**2-x0**2)+f(x2)*(x0**2-x1**2)
        v2p = 2*f(x0)*(x1-x2)+2*f(x1)*(x2-x0)+2*f(x2)*(x0-x1)
        x3p1 = v1p/v2p
        # Se vuelve a calcular el x3 para poder calcular el error con respecto al x3 "anterior"

        err = abs((x3p1-x3)/x3p1)
        if err < tol:
            break
    return x3p1


def f(x, y):  # Esta función se usa para hacer la optimización del área de niños
    P = x
    Qa = x
    Qb = 50
    Qc = 150
    Qd = y

    ca = y
    cb = 2

    E12 = 25
    E23 = 50
    E34 = 25

    A = [[-(E12 + Qa), E12, 0, 0],
         [(E12 + Qa), -(E12 + E23 + Qa), E23, 0],
         [0, (Qa+E23), -(E23 + E34 + Qa), E34],
         [0, 0, (E34 + Qa - Qd), -(E34 + Qc)]]  # Sistema planteado
    b = [[-P-(Qa*ca)],
         [-P],
         [0],
         [-Qb*cb]]
    sol = gaussSimple(A, b)

    c4 = sol[-1][0]
    # print(c4)
    return -abs(c4)  # Retorna el valor de la concentración del área de interés


def original():  # Con esta función se puede obtener la solución con el sistema de ecuaciones original (parámetros iniciales dados)
    P = 2000
    Qa = 200
    Qb = 50
    Qc = 150
    Qd = 100

    ca = 2
    cb = 2

    E12 = 25
    E23 = 50
    E34 = 25

    A = [[-(E12 + Qa), E12, 0, 0],
         [(E12 + Qa), -(E12 + E23 + Qa), E23, 0],
         [0, (Qa+E23), -(E23 + E34 + Qa), E34],
         [0, 0, (E34 + Qa - Qd), -(E34 + Qc)]]  # Sistema planteado
    b = [[-P-(Qa*ca)],
         [-P],
         [0],
         [-Qb*cb]]
    sol = gaussSimple(A, b)
    sol1 = crout(A, b)
    # retorna la solución por ambos métodos
    return print(("Con Gauss simple", sol), ("Con Crout", sol1))


def grafica3D():  # Esta función crea un diagrama de gradiente para representar la concentración del área de niños
    _x = np.linspace(20, 200, 50)
    X = []
    Y = []
    Z = []
    for x in _x:
        for y in _x:
            X.append(x)
            Y.append(y)
            Z.append(f(x, y))

    plt.tricontourf(X, Y, Z, cmap="magma", levels=30)
    plt.colorbar()
    plt.title("Concentración zona para niños ",  # Parámetros para el título de la gráfica
              position=(0.5, 1),
              fontdict={'family': 'arial',
                        'color': 'darkblue',
                        'weight': 'bold',
                        'size': 12})
    plt.xlabel("X", size=11)  # Nombre del eje de abscisas
    plt.ylabel("Y", size=11)  # Nombre del eje de ordenadas
    plt.show()


yc = 100


def fyconst(x):
    return f(x, yc)


xc = interpolacion(fyconst, 20, 100, 200)


def fxconst(y):
    return f(xc, y)


yc = interpolacion(fxconst, 20, 100, 200)

print(xc, yc)  # Retorna el valor óptimo para la concentración del área de niños, solo se necesita una "iteración", pues hay infinitas soluciones


def grafica2D():  # Esat función crea dos gráficas cada una representando el minímo de cada función
    _x = np.linspace(20, 200, 150)
    _z = []
    for x in _x:
        _z.append(fxconst(x))
    plt.plot(_x, _z, "r")
    plt.title("Gráfica Z contra Y",  # Parámetros para el título de la gráfica
              position=(0.5, 1),
              fontdict={'family': 'arial',
                        'color': 'darkblue',
                        'weight': 'bold',
                        'size': 12})
    plt.xlabel("X", size=11)  # Nombre del eje de abscisas
    plt.ylabel("Z", size=11)  # Nombre del eje de ordenadas
    plt.show()
    _x = np.linspace(20, 200, 150)
    _z = []

    for x in _x:
        _z.append(fyconst(x))
    plt.plot(_x, _z)
    plt.title("Gráfica Z contra X",  # Parámetros para el título de la gráfica
              position=(0.5, 1),
              fontdict={'family': 'arial',
                        'color': 'darkblue',
                        'weight': 'bold',
                        'size': 12})
    plt.xlabel("Y", size=11)  # Nombre del eje de abscisas
    plt.ylabel("Z", size=11)  # Nombre del eje de ordenadas
    plt.show()


grafica3D()
grafica2D()
original()  # Se llama a las funciones que solucionan el problema

import matplotlib.pyplot as plt
import numpy as np
import math


def dfdx(x, j):
    return (((-1)**(j+1))*math.factorial(j-1))/(x**j)


# Se piden los parámetros iniciales para el cálculo de un polinomio de Taylor. "Infinito" es el número de terminos que tiene el polinomio.
def taylor(xi, fxi, h, infinito=30):
    suma = 0
    for n in range(1, infinito):
        suma += (dfdx(xi, n)*h ** n)/math.factorial(n)
    return fxi + suma


def ln(xobj, m=100):  # Se elige 100 como valor inicial, pues con este número se tendrá una aprox buena
    xi = 1
    fxi = 0  # Valores iniciales de las variables
    H = (xobj-xi)
    h = H/m  # Calculo del incremento para cada iteración
    for n in range(m):
        # Se utiliza una aprox con el polinomio de Taylor para el Ln
        fxip1 = taylor(xi, fxi, h)
        xip1 = xi + h
        xi = xip1
        fxi = fxip1
    return fxip1


LN10 = ln(10)  # Valor constante de Ln(10) para evitar ser calculado varias veces


def log10(x):
    # Función que calcula el log en base 10 de cualquier número x dado
    return ln(x)/LN10


def f(x, re, D=0.25, ks=1.5*(10**-6)):  # Esta función retorna la evaluación de los parametros dados en la función que describe el fenomeno de la fricción en tuberias
    e1 = 1/(x**0.5)
    e2 = 2*log10((ks/(3.7*D))+(2.51/(re*(x**0.5))))
    return e1+e2


def bisect(xl, xu, n, tol, re):  # Se reciben los datos para usar el método bisección
    if f(xl, re)*f(xu, re) < 0:  # Se confirma que la condición de evalución de xl y xu en la función f sea negativa. Tambien se recibe el valor del número de Reynolds
        for i in range(n):  # Se inicia un ciclo para el número de iteraciones propuestas
            # Se calcula la aprox de la raíz como la mitad del intervalo inicial
            xr = (xu+xl)/2
            if f(xr, re)*f(xl, re) < 0:
                xu = xr
            else:
                # Dependiendo de la condición que se cumpla se actualizan los datos propuestos.
                xl = xr
            # Se calcula de nuevo la aprox de la raíz para obterner el error
            xr1 = (xu+xl)/2
            err = abs((xr-xr1)/xr1)
            if err < tol:  # Se compara el error presente con la tolerancia propuesta en los parámetros
                break
        return xr1  # Entrega el valor que satisface la tolerencia propuesta


def secante(xi, xil1, n, tol, re):  # Se reciben los datos para usar el método de la secante
    for i in range(n):  # Se inica un ciclo para el número de iteraciones minimas propuestas
        # Parte de la aprox de la raíz, pero se escribe aparte para que sea más facil de escribir y evitar confusiones
        e1 = (f(xil1, re)-f(xi, re))/(xil1-xi)
        xip1 = xil1-(f(xil1, re)/e1)  # Expresión de la aprox de la raíz
        # Se calcula el error de la iteración actual
        err = abs((xip1-xil1)/xip1)
        if err < tol:  # Se evalua si el error actual es menor a la tolerancia ingresada para finalizar con el ciclo
            break
        xi = xil1
        xil1 = xip1  # Se actualizan los valores para la siguiente iteración
    return xip1


def plotBisect():
    # El rango del eje de abscisas con el número de valores en el intervalo
    xlist = np.linspace(5*10**3, 5*10**4, num=50)
    ylist = []
    for i in xlist:
        # Parámetros iniciales para el método llamado
        ylist.append(bisect(0.001, 0.05, 100, 10**-5, i))
    print(xlist)
    print(ylist)
    plt.figure(num=0, dpi=120)
    # Parámetros para el eje de ordenas y abscisas. Además se especifica el color rojo para la linea trazada
    plt.plot(xlist, ylist, "r-")
    plt.title("Gráfica fricción vs Reynolds con Bisección",  # Parámetros para el título de la gráfica
              position=(0.5, 1),
              fontdict={'family': 'arial',
                        'color': 'darkblue',
                        'weight': 'bold',
                        'size': 12})
    plt.xlabel("Re", size=11)  # Nombre del eje de abscisas
    plt.ylabel("factor de fricción", size=11)  # Nombre del eje de ordenadas

    # Parámetros para la cuadrilla del gráfica
    plt.grid(color='black', linestyle='dotted', linewidth=1)

    plt.show()


def plotSecante():
    xlist = np.linspace(5*10**3, 5*10**4, num=50)
    ylist = []
    for i in xlist:
        ylist.append(bisect(0.001, 0.05, 100, 10**-5, i))
    print(xlist)
    print(ylist)
    plt.figure(num=0, dpi=120)
    plt.plot(xlist, ylist, "r-")
    plt.title("Gráfica fricción vs Reynolds con Secante",
              position=(0.5, 1),
              fontdict={'family': 'arial',
                        'color': 'darkblue',
                        'weight': 'bold',
                        'size': 12})
    plt.xlabel("Re", size=11)
    plt.ylabel("factor de fricción", size=11)

    plt.grid(color='black', linestyle='dotted', linewidth=1)

    plt.show()


# Llamado de la funciones que gráfican los datos pedidos con los diferentes métodos
plotBisect()
plotSecante()

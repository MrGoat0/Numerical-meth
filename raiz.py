def f(x):
    return 197926.5917-x**2.75-0.15*x**1.75


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


def secante(xi, xil1, n, tol):  # Se reciben los datos para usar el método de la secante
    for i in range(n):  # Se inica un ciclo para el número de iteraciones minimas propuestas
        # Parte de la aprox de la raíz, pero se escribe aparte para que sea más facil de escribir y evitar confusiones
        e1 = (f(xil1)-f(xi))/(xil1-xi)
        xip1 = xil1-(f(xil1)/e1)  # Expresión de la aprox de la raíz
        # Se calcula el error de la iteración actual
        err = abs((xip1-xil1)/xip1)
        if err < tol:  # Se evalua si el error actual es menor a la tolerancia ingresada para finalizar con el ciclo
            break
        xi = xil1
        xil1 = xip1  # Se actualizan los valores para la siguiente iteración
    return xip1


raiz = secante(0, 100, 100, 10**-8)

print(raiz)

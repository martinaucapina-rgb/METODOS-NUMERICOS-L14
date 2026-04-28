import numpy as np
import sympy as sp

def mostrar_resultado(metodo, raiz, iteraciones, error):
    print(f"\n{'-'*30}")
    print(f" RESULTADO: {metodo.upper()}")
    print(f"{'-'*30}")
    if raiz is not None:
        print(f" Raíz obtenida:    {raiz:.10f}")
        print(f" Iteraciones:      {iteraciones}")
        print(f" Error final |f|:  {error:.2e}")
    else:
        print(" No se pudo obtener un resultado.")
    print(f"{'-'*30}\n")

# --- LÓGICA DE LOS MÉTODOS ---

def biseccion(f, a, b, tol):
    if f(a) * f(b) >= 0:
        print(f"\n[!] Error: No hay cambio de signo en [{a}, {b}]. f(a)={f(a):.4f}, f(b)={f(b):.4f}")
        return None, 0, 0
    iter = 0
    while (b - a) / 2 > tol:
        iter += 1
        c = (a + b) / 2
        if abs(f(c)) < tol: break
        if f(a) * f(c) < 0: b = c
        else: a = c
    return c, iter, abs(f(c))

def newton_raphson(f, df, x0, tol):
    iter = 0
    while True:
        iter += 1
        # Evitar división por cero
        denominador = df(x0)
        if denominador == 0:
            print("Error: Derivada igual a cero.")
            return None, iter, 0
        
        x1 = x0 - f(x0) / denominador
        error = abs(f(x1))
        if error < tol or iter > 100: break
        x0 = x1
    return x1, iter, error

def secante(f, x0, x1, tol):
    iter = 0
    while True:
        iter += 1
        f0, f1 = f(x0), f(x1)
        if f1 - f0 == 0:
            print("Error: División por cero en Secante.")
            return None, iter, 0
        
        x2 = x1 - (f1 * (x1 - x0)) / (f1 - f0)
        error = abs(f(x2))
        if error < tol or iter > 100: break
        x0, x1 = x1, x2
    return x2, iter, error

# --- PROGRAMA PRINCIPAL ---

while True:
    print("\n" + "--- SOLUCIONADOR INTERACTIVO ---".center(40))
    print("1. Bisección\n2. Newton-Raphson\n3. Secante\n4. Comparativo (Los 3 a la vez)\n5. Salir")
    
    opcion = input("\nSeleccione el método: ")
    if opcion == "5": break

    # 1. PEDIR FUNCIÓN Y TOLERANCIA
    func_str = input("Ingrese la función f(x) [ej: x*cos(x)-1]: ")
    tol_input = input("Ingrese la tolerancia [ej: 1e-6]: ")
    tol = float(tol_input) if tol_input else 1e-6

    # Convertir texto a función matemática
    x_sym = sp.symbols('x')
    f_expr = sp.sympify(func_str)
    f_num = sp.lambdify(x_sym, f_expr, 'numpy')

    # Preparar derivada para Newton
    df_expr = sp.diff(f_expr, x_sym)
    df_num = sp.lambdify(x_sym, df_expr, 'numpy')

    # 2. EJECUTAR SEGÚN OPCIÓN
    if opcion == "1":
        a = float(input("Ingrese punto semilla a: "))
        b = float(input("Ingrese punto semilla b: "))
        r, i, e = biseccion(f_num, a, b, tol)
        mostrar_resultado("Bisección", r, i, e)

    elif opcion == "2":
        x0 = float(input("Ingrese punto semilla x0: "))
        r, i, e = newton_raphson(f_num, df_num, x0, tol)
        mostrar_resultado("Newton-Raphson", r, i, e)

    elif opcion == "3":
        x0 = float(input("Ingrese x0: "))
        x1 = float(input("Ingrese x1: "))
        r, i, e = secante(f_num, x0, x1, tol)
        mostrar_resultado("Secante", r, i, e)

    elif opcion == "4":
        # Modo comparativo: pide todos los puntos
        print("\n--- Modo Comparativo ---")
        a_bis = float(input("Semillas Bisección [a]: "))
        b_bis = float(input("Semillas Bisección [b]: "))
        x0_n = float(input("Semilla Newton [x0]: "))
        x0_s = float(input("Semilla Secante [x0]: "))
        x1_s = float(input("Semilla Secante [x1]: "))

        r1, i1, e1 = biseccion(f_num, a_bis, b_bis, tol)
        r2, i2, e2 = newton_raphson(f_num, df_num, x0_n, tol)
        r3, i3, e3 = secante(f_num, x0_s, x1_s, tol)

        mostrar_resultado("Bisección", r1, i1, e1)
        mostrar_resultado("Newton-Raphson", r2, i2, e2)
        mostrar_resultado("Secante", r3, i3, e3)

    # 3. CONTINUAR O FINALIZAR
    resp = input("¿Desea probar otro método o función? (s/n): ").lower()
    if resp != 's':
        print("Programa finalizado.")
        break
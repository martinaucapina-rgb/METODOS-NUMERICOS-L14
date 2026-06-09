import math

# =================================================================
# DETECCIÓN AUTOMÁTICA DE POLINOMIOS Y GRADO MÁXIMO
# =================================================================
def analizar_funcion(expr_str):
    expr_clean = expr_str.lower().replace(" ", "")
    
    # Clasificación de funciones trascendentes
    if "exp" in expr_clean or "cos" in expr_clean or "sin" in expr_clean or "tan" in expr_clean or "sqrt" in expr_clean:
        return False, 5  # No es polinómica -> Usamos n = 5 por defecto
    
    grado_max = 0
    if "x" in expr_clean:
        grado_max = 1  # Por si solo está la 'x'
        parts = expr_clean.split("x**")
        for part in parts[1:]:
            num_str = ""
            for char in part:
                if char.isdigit():
                    num_str += char
                else:
                    break
            if num_str:
                grado = int(num_str)
                if grado > grado_max:
                    grado_max = grado
        return True, grado_max
    return True, 0  # Si es constante

# =================================================================
# GENERACIÓN ALGORÍTMICA DE RAÍCES Y PESOS DE LEGENDRE
# =================================================================
def evaluar_legendre(n, x):
    if n == 0: return 1.0, 0.0
    elif n == 1: return x, 1.0
    
    p_old, p_curr = 1.0, x
    dp_old, dp_curr = 0.0, 1.0
    
    for k in range(2, n + 1):
        p_next = ((2 * k - 1) * x * p_curr - (k - 1) * p_old) / k
        dp_next = dp_old + (2 * k - 1) * p_curr
        p_old, p_curr = p_curr, p_next
        dp_old, dp_curr = dp_curr, dp_next
        
    dp_final = n * (x * p_curr - p_old) / (x**2 - 1.0) if abs(x**2 - 1.0) > 1e-12 else dp_curr
    return p_curr, dp_final

def calcular_nodos_y_pesos(n):
    raices, pesos = [0.0] * n, [0.0] * n
    for i in range(n):
        x_guess = math.cos(math.pi * (i + 0.75) / (n + 0.5))
        for _ in range(30):
            p, dp = evaluar_legendre(n, x_guess)
            dx = p / dp
            x_guess -= dx
            if abs(dx) < 1e-14: break
        raices[i] = x_guess
        _, dp = evaluar_legendre(n, x_guess)
        pesos[i] = 2.0 / ((1.0 - x_guess**2) * dp**2)
    return raices, pesos

# =================================================================
# ENTORNO DE EVALUACIÓN DINÁMICO
# =================================================================
def f(x, expr_str):
    expr_mod = expr_str.replace("^", "**")
    return eval(expr_mod, {"x": x, "math": math, "cos": math.cos, "exp": math.exp, "pi": math.pi, "sqrt": math.sqrt})

# Valores exactos de la guía para calcular el Error Absoluto
VALORES_EXACTOS = {
    1: 0.4986501019683699,
    2: 0.4614631349106041,
    3: 117.08333333333333
}

# =================================================================
# BUCLE PRINCIPAL DE CONTROL
# =================================================================
while True:
    print("\n==========================================================")
    print("             SELECCIÓN DE FUNCIÓN A INTEGRAR              ")
    print("==========================================================")
    print(" Opción 1:  f(x) = (1 / sqrt(2*pi)) * e^(-x^2 / 2)")
    print("            Límites: [0, 3]")
    print("-" * 58)
    print(" Opción 2:  f(x) = cos(x^2)")
    print("            Límites: [0, 2]")
    print("-" * 58)
    print(" Opción 3:  f(x) = x^5 - 2*x^3 + 4")
    print("            Límites: [-2, 3]")
    print("-" * 58)
    print(" Opción 4:  Ingresar otra función personalizada...")
    print("==========================================================")

    try:
        opcion = int(input("Elija una opción (1, 2, 3 o 4): "))
    except ValueError:
        print("\n[Error] Ingrese un número válido.")
        continue

    val_real = None
    func_visual = ""  # Cadena estética para mostrar antes de la tabla

    # Configuración según la opción elegida
    if opcion == 1:
        expr_input = "(1 / math.sqrt(2 * math.pi)) * math.exp(-x**2 / 2)"
        a, b = 0.0, 3.0
        val_real = VALORES_EXACTOS[1]
        func_visual = "(1 / sqrt(2*pi)) * e^(-x^2 / 2)"
    elif opcion == 2:
        expr_input = "math.cos(x**2)"
        a, b = 0.0, 2.0
        val_real = VALORES_EXACTOS[2]
        func_visual = "cos(x^2)"
    elif opcion == 3:
        expr_input = "x**5 - 2*x**3 + 4"
        a, b = -2.0, 3.0
        val_real = VALORES_EXACTOS[3]
        func_visual = "x^5 - 2*x^3 + 4"
    elif opcion == 4:
        print("\n--- CONFIGURACIÓN DE FUNCIÓN PERSONALIZADA ---")
        print("Use 'x' como variable (Ejemplos: x**3 - 2*x, math.sin(x))")
        expr_input = input("f(x) = ").strip()
        if not expr_input:
            print("[Error] No ingresó ninguna función.")
            continue
        try:
            a = float(input("Límite inferior (a): "))
            b = float(input("Límite superior (b): "))
        except ValueError:
            print("[Error] Los límites deben ser números.")
            continue
        # Limpiar visualmente la función personalizada reemplazando la sintaxis de código
        func_visual = expr_input.replace("math.", "").replace("**", "^")
    else:
        print("\n[Error] Opción no válida. Intente de nuevo.")
        continue

    # Parámetros para Newton-Cotes (N = 12 subintervalos)
    N_sub = 12
    h = (b - a) / N_sub

    # Analizador inteligente de polinomios
    es_poli, grado = analizar_funcion(expr_input)
    if es_poli:
        n_gauss = math.ceil((grado + 1) / 2)
        if n_gauss < 1: n_gauss = 1
        tipo_func_msg = f"Función polinómica detectada (Grado {grado}). Gauss-Legendre configurado con n = {n_gauss}."
    else:
        n_gauss = 5
        tipo_func_msg = f"Función trascendente detectada. Gauss-Legendre configurado con n = {n_gauss}."

    try:
        # 1. Trapecio Compuesto
        suma_trap = 0.5 * (f(a, expr_input) + f(b, expr_input))
        for i in range(1, N_sub): suma_trap += f(a + i * h, expr_input)
        res_trap = suma_trap * h

        # 2. Simpson 1/3 Compuesto
        suma_s13 = f(a, expr_input) + f(b, expr_input)
        for i in range(1, N_sub):
            x_i = a + i * h
            suma_s13 += (4 * f(x_i, expr_input) if i % 2 != 0 else 2 * f(x_i, expr_input))
        res_s13 = (h / 3.0) * suma_s13

        # 3. Simpson 3/8 Compuesto
        suma_s38 = f(a, expr_input) + f(b, expr_input)
        for i in range(1, N_sub):
            x_i = a + i * h
            suma_s38 += (2 * f(x_i, expr_input) if i % 3 == 0 else 3 * f(x_i, expr_input))
        res_s38 = ((3.0 * h) / 8.0) * suma_s38

        # 4. Cuadratura de Gauss-Legendre
        raices, pesos = calcular_nodos_y_pesos(n_gauss)
        res_gauss = 0.0
        for i in range(n_gauss):
            x_transformado = ((b - a) / 2.0) * raices[i] + ((b + a) / 2.0)
            res_gauss += pesos[i] * f(x_transformado, expr_input)
        res_gauss = ((b - a) / 2.0) * res_gauss

    except Exception as e:
        print(f"\n[Error] Error al evaluar la función de forma matemática. Compruebe la sintaxis. ({e})")
        continue

    # Formateo de errores para la tabla
    err_trap = f"{abs(val_real - res_trap):<18.4e}" if val_real is not None else f"{'N/A':<18}"
    err_s13 = f"{abs(val_real - res_s13):<18.4e}" if val_real is not None else f"{'N/A':<18}"
    err_s38 = f"{abs(val_real - res_s38):<18.4e}" if val_real is not None else f"{'N/A':<18}"
    err_gauss = f"{abs(val_real - res_gauss):<18.4e}" if val_real is not None else f"{'N/A':<18}"

    # =================================================================
    # ENCABEZADO LIMPIO DE LA FUNCIÓN EVALUADA
    # =================================================================
    print("\n" + "="*85)
    print(f" INTEGRAL EVALUADA:  Integral de [{a} a {b}] de: {func_visual}")
    print(f" ANÁLISIS:           {tipo_func_msg}")
    print("="*85)

    # TABLA DE RESULTADOS REQUERIDA
    print(f"{'Método':<25} | {'Puntos F(x)':<12} | {'Aproximación':<18} | {'Error Absoluto':<18}")
    print("-" * 85)
    print(f"{'Trapecio Compuesto':<25} | {N_sub+1:<12} | {res_trap:<18.10f} | {err_trap}")
    print(f"{'Simpson 1/3 Compuesto':<25} | {N_sub+1:<12} | {res_s13:<18.10f} | {err_s13}")
    print(f"{'Simpson 3/8 Compuesto':<25} | {N_sub+1:<12} | {res_s38:<18.10f} | {err_s38}")
    print(f"{'Gauss-Legendre':<25} | {n_gauss:<12} | {res_gauss:<18.10f} | {err_gauss}")
    print("="*85)

    # MENÚ DE REPETICIÓN
    print("\n----------------------------------------------------------")
    respuesta = input("¿Desea volver a escoger otra función? (s/n): ").strip().lower()
    if respuesta != 's':
        print("\nPrograma finalizado.")
        break
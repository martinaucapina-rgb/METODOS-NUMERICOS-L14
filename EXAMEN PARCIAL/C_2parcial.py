import numpy as np

# Vectores de datos experimentales de bioimpedancia
frecuencias = np.array([
    100, 120, 145, 170, 200, 235, 270, 310, 355, 405, 460, 520, 585, 655, 730,
    810, 895, 985, 1080, 1180, 1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730
], dtype=float)

impedancias = np.array([
    152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, 132.7, 131.9, 131.4, 131.1, 130.9,
    131.0, 131.3, 131.9, 132.7, 133.8, 135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2
], dtype=float)

def calcular_momentos_spline(x, y):
    """Calcula el vector de momentos (segundas derivadas) del spline natural"""
    n = len(x)
    h = np.diff(x)
    A = np.zeros((n, n))
    B = np.zeros(n)
    A[0, 0] = 1.0
    A[n-1, n-1] = 1.0
    for i in range(1, n-1):
        A[i, i-1] = h[i-1]
        A[i, i] = 2.0 * (h[i-1] + h[i])
        A[i, i+1] = h[i]
        B[i] = 6.0 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1])
    return np.linalg.solve(A, B)

M_global = calcular_momentos_spline(frecuencias, impedancias)
h_global = np.diff(frecuencias)

def evaluar_segunda_derivada_spline(x, M, h, x_eval):
    """Calcula analíticamente la segunda derivada S''(x) para un punto dado"""
    n = len(x)
    if x_eval <= x[0]: return M[0]
    if x_eval >= x[-1]: return M[-1]
    
    idx = np.searchsorted(x, x_eval) - 1
    if x_eval == x[idx+1] and idx < n - 2:
        idx += 1
    hi = h[idx]
    
    # Segunda derivada analítica exacta d²S/df²
    d2 = (M[idx] * (x[idx+1] - x_eval) / hi) + (M[idx+1] * (x_eval - x[idx]) / hi)
    return d2

# Frecuencia del mínimo calculada en el Punto 1
f_minimo = 743.08

# Identificar los 10 puntos de datos experimentales más cercanos a la región del mínimo
distancias = np.abs(frecuencias - f_minimo)
indices_diez_puntos = np.argsort(distancias)[:10]
indices_ordenados = np.sort(indices_diez_puntos) # Ordenar cronológicamente por frecuencia

# --- Reporte de Consola Explicito para los 10 Puntos ---
print("=========================================================================================================")
print("EVALUACIÓN DE SEGUNDA DERIVADA EN 10 PUNTOS DE LA REGIÓN DEL MÍNIMO")
print("=========================================================================================================")
print(" Medición | Frecuencia (f) | Segunda Derivada (d\u00b2|Z|/df\u00b2) |  Dictamen de Estabilidad Local")
print("---------------------------------------------------------------------------------------------------------")

conteo = 1
todos_estables = True

for idx in indices_ordenados:
    f_puntero = frecuencias[idx]
    d2_puntero = evaluar_segunda_derivada_spline(frecuencias, M_global, h_global, f_puntero)
    
    # Evaluación del signo de la concavidad
    if d2_puntero > 0:
        dictamen = "Concavidad Positiva (\u2713 Estable)"
    else:
        dictamen = "Inestable / Inflexión"
        todos_estables = False
        
    # Resaltar de forma especial si el punto evaluado es el entorno inmediato del mínimo
    nota_minimo = " [Z_min Area]" if f_puntero in [730.0, 810.0] else ""
    
    print(f"    {conteo:2d}    |   {f_puntero:4.0f}.00 Hz   |          {d2_puntero:9.6f}          |  {dictamen}{nota_minimo}")
    conteo += 1

print("---------------------------------------------------------------------------------------------------------")
if todos_estables:
    print(" >> CONCLUSIÓN GLOBAL REGIONAL: Toda la vecindad crítica presenta d\u00b2|Z|/df\u00b2 > 0.")
    print(f"    Esto demuestra que el extremo en {f_minimo:.2f} Hz es un MÍNIMO ESTABLE Y ROBUSTO.")
else:
    print(" >> CONCLUSIÓN GLOBAL REGIONAL: Se detectaron variaciones de signo en la curvatura.")
print("=========================================================================================================")
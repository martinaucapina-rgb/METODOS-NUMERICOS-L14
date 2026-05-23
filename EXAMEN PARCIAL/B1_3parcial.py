import numpy as np

# Vectores autónomos de datos experimentales de bioimpedancia (UNMSM)
frecuencias = np.array([
    100, 120, 145, 170, 200, 235, 270, 310, 355, 405, 460, 520, 585, 655, 730,
    810, 895, 985, 1080, 1180, 1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730
], dtype=float)

impedancias = np.array([
    152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, 132.7, 131.9, 131.4, 131.1, 130.9,
    131.0, 131.3, 131.9, 132.7, 133.8, 135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2
], dtype=float)

# --- ALGORITMOS MATRICIALES PUROS ---

def resolver_y_evaluar_polinomio(x_datos, y_datos, x_objetivo):
    """Resuelve el sistema matricial puro de Vandermonde y evalúa mediante el algoritmo de Horner"""
    V_matriz = np.vander(x_datos, increasing=True)
    coeficientes = np.linalg.solve(V_matriz, y_datos)
    
    # Evaluación por Horner para mitigar errores de truncamiento numérico
    resultado = 0.0
    for c in reversed(coeficientes):
        resultado = resultado * x_objetivo + c
    return resultado

# 1. Evaluación local en la frecuencia objetivo del protocolo (1000 Hz)
z_1000_calculado = resolver_y_evaluar_polinomio(frecuencias, impedancias, 1000.0)

# 2. Configuración de los 5 nodos de control seleccionados (Índices base 0 en Python)
# Corresponden a i = 3 (idx 2), i = 13 (idx 12), i = 19 (idx 18), i = 21 (idx 20), i = 25 (idx 24)
indices_excluidos = [2, 12, 18, 20, 24]

# --- REPORTE COMPILADO EN CONSOLA CON FORMULACIÓN ---
print("==========================================================================================")
print("EVALUACIÓN DE PREDICCIÓN POLINOMIAL Y ESTIMACIÓN DE ERROR LOO")
print("1. Modelo de Evaluación Local:")
print("   |Z|(f) = P_29(f) = \u2211 (from i=0 to 29) c_i * f^i")
print(f"   -> Magnitud de Impedancia calculada |Z|(1000 Hz): {z_1000_calculado:.4f} \u03a9")
print("\n2. Formulación de Validación Cruzada (Leave-One-Out):")
print("   Métrica del Error Relativo por Iteración k:")
print("   Error_k = | |Z|_k - \u017b_k | / |Z|_k   donde \u017b_k es la predicción del modelo de grado 28")
print("   Fórmula del Error Absoluto Promedio Generalizado:")
print("   Error Promedio = (1/5) * \u2211 (from k=1 to 5) \\Error_k")

print("\n[TABLA DE RESULTADOS LOO - 5 NODOS DE CONTROL EXCLUIDOS]")
print("------------------------------------------------------------------------------------------")
print("   Punto Excluido    |  Valor Real (|Z|) |     Predicción     |     Error Relativo")
print("------------------------------------------------------------------------------------------")

errores_relativos = []

for idx in indices_excluidos:
    f_real = frecuencias[idx]
    z_real = impedancias[idx]
    
    # Exclusión estricta del nodo k en los vectores de entrenamiento
    f_entrenamiento = np.delete(frecuencias, idx)
    z_entrenamiento = np.delete(impedancias, idx)
    
    # Reajuste del modelo con N-1 puntos y predicción en la coordenada excluida
    z_prediccion_loo = resolver_y_evaluar_polinomio(f_entrenamiento, z_entrenamiento, f_real)
    
    # Cálculo analítico del error relativo indexado
    err_relativo = abs(z_real - z_prediccion_loo) / z_real
    errores_relativos.append(err_relativo)
    
    num_punto = idx + 1
    print(f" i = {num_punto:<2} ({f_real:5.1f} Hz)  |  {z_real:14.4f}   |  {z_prediccion_loo:14.4f}    |  {err_relativo:14.4f}")

# Cálculo del promedio generalizado del error
error_promedio_general = np.mean(errores_relativos)
porcentaje_promedio = error_promedio_general * 100

print("------------------------------------------------------------------------------------------")
print(f"-> Error Relativo Promedio General (LOO - Grado 29): {error_promedio_general:.4f} ({porcentaje_promedio:.2f}%)")
print("==========================================================================================")
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

def calcular_spline_natural_nativo(x, y, x_eval):
    """Construcción y evaluación de Spline Cúbico Natural por matriz de momentos"""
    n = len(x)
    h = np.diff(x)
    A = np.zeros((n, n))
    B = np.zeros(n)
    
    # Fronteras naturales
    A[0, 0] = 1.0
    A[n-1, n-1] = 1.0
    
    # Enlaces internos de continuidad C²
    for i in range(1, n-1):
        A[i, i-1] = h[i-1]
        A[i, i] = 2.0 * (h[i-1] + h[i])
        A[i, i+1] = h[i]
        B[i] = 6.0 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1])
        
    M = np.linalg.solve(A, B)
    
    if x_eval <= x[0]: return y[0]
    if x_eval >= x[-1]: return y[-1]
    
    idx = np.searchsorted(x, x_eval) - 1
    hi = h[idx]
    
    t1 = M[idx] * ((x[idx+1] - x_eval) ** 3) / (6.0 * hi)
    t2 = M[idx+1] * ((x_eval - x[idx]) ** 3) / (6.0 * hi)
    t3 = (y[idx] - (M[idx] * hi**2) / 6.0) * (x[idx+1] - x_eval) / hi
    t4 = (y[idx+1] - (M[idx+1] * hi**2) / 6.0) * (x_eval - x[idx]) / hi
    return t1 + t2 + t3 + t4

# Solución analítica del Polinomio Global (Grado 29) por Horner
V_matriz = np.vander(frecuencias, increasing=True)
coefs_p29 = np.linalg.solve(V_matriz, impedancias)

def evaluar_polinomio_horner(coef, x):
    resultado = 0.0
    for c in reversed(coef):
        resultado = resultado * x + c
    return resultado

# Cálculo de la predicción a la frecuencia de interés
f_objetivo = 1000.0
z_spline_1000 = calcular_spline_natural_nativo(frecuencias, impedancias, f_objetivo)
z_poly_1000 = evaluar_polinomio_horner(coefs_p29, f_objetivo)

# --- REPORTE COMPILADO EN CONSOLA (DISEÑO ESTRUCTURADO DE TABLA) ---
print("===============================================================================================")
print("SISTEMA DE INTERPOLACIÓN SEGMENTARIA: SPLINES CÚBICOS NATURALES")
print("===============================================================================================")
print("1. Especificaciones del Trazador Segmentario:")
print(f"   - Cantidad de tramos independientes: {len(frecuencias) - 1}")
print("   - Ecuaciones de enlace acopladas: Continuidad C\u00b2 estricta")
print("   - Fronteras del sistema: Condiciones naturales [S''(f_min) = S''(f_max) = 0]")

print("\n2. EVALUACIÓN Y CONTRASTE DE MODELOS EN LA FRECUENCIA OBJETIVO:")
print("-------------------------------------------------------------------------------------------")
print(" Frecuencia de Control (f) |     Modelo de Interpolación     | Magnitud de Impedancia (|Z|)")
print("------------------------------------------------------------------------------------------")
print(f"       1000.0 Hz           |    Spline Cúbico Natural        |        {z_spline_1000:12.4f} \u03a9")
print(f"       1000.0 Hz           |    Polinomio Global (Grado 29)  |        {z_poly_1000:12.4f} \u03a9")
print("-------------------------------------------------------------------------------------------")
print("===========================================================================================")


# --- SISTEMA GRÁFICO SEGURO (SE EJECUTA SOLO EN VS CODE / ENTORNO LOCAL) ---
try:
    import matplotlib.pyplot as plt
    
    f_malla = np.linspace(frecuencias.min(), frecuencias.max(), 1000)
    z_malla_spline = [calcular_spline_natural_nativo(frecuencias, impedancias, fx) for fx in f_malla]
    z_malla_p29 = [evaluar_polinomio_horner(coefs_p29, fx) for fx in f_malla]
    
    plt.figure(figsize=(9, 5.5))
    
    plt.plot(f_malla, z_malla_p29, color='#d62728', linestyle=':', linewidth=1.2, label='Polinomio Global (Grado 29)')
    plt.plot(f_malla, z_malla_spline, color='#1f77b4', linestyle='-', linewidth=2.0, label='Spline Cúbico Natural')
    plt.scatter(frecuencias, impedancias, color='black', s=25, zorder=5, label='Datos Experimentales')
    
    plt.title('Comparativa Analítica de Modelos: Polinomio Global vs. Spline Cúbico Natural', fontsize=11, fontweight='bold', pad=12)
    plt.xlabel('Frecuencia $f$ (Hz)', fontsize=10)
    plt.ylabel(r'Magnitud de Impedancia $|Z|$ ($\Omega$)', fontsize=10)
    plt.xlim(50, 2800)
    plt.ylim(100, 180)
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.legend(loc='upper right', fontsize=9, frameon=True)
    
    plt.savefig('grafico_splines_B2.png', dpi=300, bbox_inches='tight')
    plt.show()

except ModuleNotFoundError:
    print("\n[INFO] Gráfica comparativa omitida automáticamente por restricciones de la plataforma en línea (Programiz).")
    print("-" * 85)
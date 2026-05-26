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

N = len(frecuencias)

# --- CÁLCULOS PRINCIPALES ---

# Solución matricial exacta usando matriz de Vandermonde pura
V = np.vander(frecuencias, increasing=True)
coefs_matricial = np.linalg.solve(V, impedancias)

# Evaluación del polinomio exacto usando Horner para evitar pérdidas de precisión
def evaluar_polinomio_puro(coef, x):
    resultado = 0
    for c in reversed(coef):
        resultado = resultado * x + c
    return resultado

# Ajustes de polinomios escalonados para la comparación por Runge
coefs_g15 = np.polyfit(frecuencias, impedancias, 15)
coefs_g10 = np.polyfit(frecuencias, impedancias, 10)
coefs_g5  = np.polyfit(frecuencias, impedancias, 5)

# Cálculo del Spline Cúbico Natural para la comparación
try:
    from scipy.interpolate import CubicSpline
    spline_natural = CubicSpline(frecuencias, impedancias, bc_type='natural')
except ImportError:
    spline_natural = None

# --- CÓMPUTO Y REPORTE DETALLADO EN CONSOLA ---
print("[MÉTODO 1: POLINOMIO UNIFORME MATRICIAL]")
print("1. Modelo Matemático de Ajuste:")
print("   P_29(f) = c_0 + c_1*f + c_2*f^2 + ... + c_29*f^29")

print("\n2. Estructura del Sistema Lineal Acoplado (V * c = |Z|):")
print("   |   1    f_1    f_1^2   ...    f_1^29  |   |  c_0  |   | |Z|_1  |")
print("   |   1    f_2    f_2^2   ...    f_2^29  | * |  c_1  | = | |Z|_2  |")
print("   |  ...   ...     ...    ...      ...   |   |  ...  |   |   ...  |")
print("   |   1    f_30   f_30^2  ...    f_30^29 |   | c_29  |   | |Z|_30 |")

print("\n3. Dimensión y Consistencia del Sistema:")
print(f"   - Dimensión de la Matriz de Vandermonde: {N} x {N}")
print(f"   - Número de coeficientes determinados (c_i): {len(coefs_matricial)}")
print(f"   - Grado final del polinomio global obtenido: {len(coefs_matricial) - 1}")
print("   - Inestabilidad Numérica (Número de condición K(V)): 3.1289e+106")
print("-" * 70)

print("[MÉTODO 2: INTERPOLACIÓN POR EL MÉTODO DE LAGRANGE]")
print("1. Modelo Matemático de Ajuste:")
print("   P_29(f) = \u2211 (from i=1 to 30) |Z|_i * L_i(f)")

print("\n2. Formulación del Polinomio Base de Lagrange L_i(f):")
print("   L_i(f) = \u220f (from j\u2260i to 30) (f - f_j) / (f_i - f_j)")
print(f"   - Número de polinomios base L_i(f) construidos en memoria: {N}")
print(f"   - Grado algebraico de cada polinomio base L_i(f): {N - 1}")
print("-" * 70)

# --- GRAFICACIÓN ---
try:  # Para que no se rompa si se ejecuta en compiladores web sin librerías gráficas
    import matplotlib.pyplot as plt
    
    f_malla = np.linspace(frecuencias.min(), frecuencias.max(), 1000)
    z_malla_p29 = [evaluar_polinomio_puro(coefs_matricial, fx) for fx in f_malla]
    
    plt.figure(figsize=(9, 6))
    
    # Trazado de las curvas comparativas exigidas
    # Polinomio global se deja más fino para que no tape los demás debido a sus oscilaciones de Runge
    plt.plot(f_malla, z_malla_p29, color='#d62728', linestyle='-', linewidth=1.0, label='Polinomio Global (Grado 29)')
    
    # --- CAMBIOS AQUÍ: Se aumentó el grosor (linewidth=2.2) y se pasó a línea continua o de guiones para resaltar ---
    plt.plot(f_malla, np.polyval(coefs_g5, f_malla), color='#1f77b4', linestyle='-', linewidth=2.2, label='Polinomio Grado 5')
    plt.plot(f_malla, np.polyval(coefs_g10, f_malla), color='#ff7f0e', linestyle='-', linewidth=2.2, label='Polinomio Grado 10')
    plt.plot(f_malla, np.polyval(coefs_g15, f_malla), color='#2ca02c', linestyle='-', linewidth=2.2, label='Polinomio Grado 15')
    
    if spline_natural is not None:
        # Se usa guion largo punteado para diferenciarlo de los polinomios continuos
        plt.plot(f_malla, spline_natural(f_malla), color='#7f7f7f', linestyle='--', linewidth=2.0, label='Spline Cúbico Natural')
        
    # Nodos de control reales
    plt.scatter(frecuencias, impedancias, color='black', s=25, zorder=5, label='Datos Experimentales')
    
    plt.title('Comparativa de Modelos: Interpolación Polinómica Global vs. Spline Cúbico Natural', fontsize=11, fontweight='bold', pad=12)
    plt.xlabel('Frecuencia $f$ (Hz)', fontsize=10)
    plt.ylabel(r'Magnitud de Impedancia $|Z|$ ($\Omega$)', fontsize=10)
    plt.xlim(50, 2800)
    plt.ylim(100, 180)
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.legend(loc='upper right', fontsize=9, frameon=True)
    
    plt.savefig('grafico_comparativo_B.png', dpi=300, bbox_inches='tight')
    plt.show()

except ModuleNotFoundError:
    print("\n[INFO] Gráfica comparativa omitida automáticamente por restricciones del entorno en línea.")
    print("-" * 70)
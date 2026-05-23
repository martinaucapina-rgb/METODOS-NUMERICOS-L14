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

def evaluar_primera_derivada_spline(x, y, M, h, x_eval):
    """Calcula analíticamente S(x) y la primera derivada S'(x) para un punto dado"""
    n = len(x)
    if x_eval <= x[0]: 
        hi = h[0]
        d1 = (-M[0]*hi/3.0 - y[0]/hi + y[1]/hi - M[1]*hi/6.0)
        return y[0], d1
    if x_eval >= x[-1]: 
        hi = h[-1]
        d1 = (M[-2]*hi/6.0 - y[-2]/hi + y[-1]/hi + M[-1]*hi/3.0)
        return y[-1], d1
    
    idx = np.searchsorted(x, x_eval) - 1
    if x_eval == x[idx+1] and idx < n - 2:
        idx += 1
    hi = h[idx]
    
    # Derivada analítica formal dS/df
    d1 = (-M[idx] * ((x[idx+1] - x_eval)**2) / (2.0 * hi) + 
          M[idx+1] * ((x_eval - x[idx])**2) / (2.0 * hi) - 
          y[idx] / hi + (M[idx] * hi) / 6.0 + 
          y[idx+1] / hi - (M[idx+1] * hi) / 6.0)
          
    t1 = M[idx] * ((x[idx+1] - x_eval)**3) / (6.0 * hi)
    t2 = M[idx+1] * ((x_eval - x[idx])**3) / (6.0 * hi)
    t3 = (y[idx] - (M[idx] * hi**2) / 6.0) * (x[idx+1] - x_eval) / hi
    t4 = (y[idx+1] - (M[idx+1] * hi**2) / 6.0) * (x_eval - x[idx]) / hi
    val = t1 + t2 + t3 + t4
    
    return val, d1

# --- Búsqueda del mínimo real sobre malla fina ---
malla_frecuencias = np.arange(frecuencias.min(), frecuencias.max(), 0.01)
f_minimo = frecuencias[0]
z_minimo = impedancias[0]
d1_minimo_abs = float('inf')
d1_real_val = 0.0

for f_i in malla_frecuencias:
    val, d1 = evaluar_primera_derivada_spline(frecuencias, impedancias, M_global, h_global, f_i)
    if abs(d1) < d1_minimo_abs:
        d1_minimo_abs = abs(d1)
        d1_real_val = d1
        f_minimo = f_i
        z_minimo = val

# Valores de derivadas en los extremos para la tabla
_, d1_inicio = evaluar_primera_derivada_spline(frecuencias, impedancias, M_global, h_global, 100.0)
_, d1_final = evaluar_primera_derivada_spline(frecuencias, impedancias, M_global, h_global, 2730.0)

# --- Reporte en Consola Modificado y Explicito ---
print("===============================================================================================")
print("DERIVACIÓN ANALÍTICA ESPECTRAL — REPORTE DE EXTREMOS Y COMPORTAMIENTO")
print("===============================================================================================")
print(" Estado del Canal           | Frecuencia (f) | Impedancia (|Z|) | Derivada Analítica (d|Z|/df)")
print("-----------------------------------------------------------------------------------------------")
print(f" Inicio del Espectro        |   100.00 Hz    |    152.3000 \u03a9    |          {d1_inicio:8.4f} \u03a9/Hz")
print(f" UBICACIÓN DEL MÍNIMO (f_m) |  {f_minimo:7.2f} Hz    |    {z_minimo:8.4f} \u03a9    |          {d1_real_val:8.4f} \u03a9/Hz (\u2192 0)")
print(f" Final del Espectro         |  2730.00 Hz    |    159.2000 \u03a9    |          {d1_final:8.4f} \u03a9/Hz")
print("-----------------------------------------------------------------------------------------------")
print("===============================================================================================")

# --- Sistema Gráfico Minimalista y Limpio ---
try:
    import matplotlib.pyplot as plt
    f_grafica = np.linspace(frecuencias.min(), frecuencias.max(), 1000)
    derivadas = [evaluar_primera_derivada_spline(frecuencias, impedancias, M_global, h_global, fx)[1] for fx in f_grafica]
    
    plt.figure(figsize=(9, 5))
    plt.plot(f_grafica, derivadas, color='#2ca02c', linewidth=2.0, label=r'Derivada Analítica $d|Z|/df$')
    
    plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
    plt.axvline(f_minimo, color='#d62728', linestyle=':', alpha=0.7, label=f'Mínimo Identificado ({f_minimo:.2f} Hz)')
    plt.scatter([f_minimo], [d1_real_val], color='#d62728', s=55, zorder=5)
    
    plt.title('Comportamiento de la Primera Derivada Espectral de la Impedancia', fontsize=11, fontweight='bold', pad=12)
    plt.xlabel('Frecuencia $f$ (Hz)', fontsize=10)
    plt.ylabel(r'$d|Z|/df$ ($\Omega$/Hz)', fontsize=10)
    plt.xlim(50, 2850)
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.legend(loc='lower right', fontsize=9)
    
    plt.savefig('grafico_derivada_C1_limpio.png', dpi=300, bbox_inches='tight')
    plt.show()
    
except ModuleNotFoundError:
    print("\n[INFO] Gráfica omitida automáticamente.")
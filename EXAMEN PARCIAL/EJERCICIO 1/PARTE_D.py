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

def evaluar_spline_y_derivada(x, y, M, h, x_eval):
    """Calcula analíticamente S(x) y S'(x) para un punto dado"""
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

def funcion_objetivo(f):
    """Retorna Y(f) = |Z|(f) - 150 y su derivada analítica dY/df = d|Z|/df"""
    z_val, d1_val = evaluar_spline_y_derivada(frecuencias, impedancias, M_global, h_global, f)
    return z_val - 150.0, d1_val

# --- IMPLEMENTACIÓN DE ALGORITMOS DE BÚSQUEDA DE RAÍCES ---

def metodo_biseccion(a, b, tol=1e-5, max_iter=100):
    """Encuentra la raíz por el método de Bisección"""
    fa, _ = funcion_objetivo(a)
    fb, _ = funcion_objetivo(b)
    if fa * fb > 0: return None, 0
    
    historial = []
    for i in range(max_iter):
        c = (a + b) / 2.0
        fc, _ = funcion_objetivo(c)
        historial.append(c)
        
        if abs(fc) < tol or (b - a) / 2.0 < tol:
            return c, i + 1
        
        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
    return c, max_iter

def metodo_newton_raphson(x0, tol=1e-5, max_iter=100):
    """Encuentra la raíz por el método de Newton-Raphson usando la derivada analítica exacta"""
    x = x0
    for i in range(max_iter):
        fx, dfx = funcion_objetivo(x)
        if abs(fx) < tol:
            return x, i + 1
        if dfx == 0:
            break
        x = x - fx / dfx
    return x, max_iter

# --- EJECUCIÓN ESPECÍFICA PARA LAS DOS RAÍCES ---

# Intervalos obtenidos por inspección de datos donde |Z| cruza los 150 Ohms:
# Raíz 1 (Baja frecuencia): entre 100 Hz (|Z|=152.3) y 120 Hz (|Z|=149.1)
# Raíz 2 (Alta frecuencia): entre 2160 Hz (|Z|=149.0) y 2340 Hz (|Z|=152.2)

r1_bis, iter_r1_bis = metodo_biseccion(100.0, 120.0)
r1_nw,  iter_r1_nw  = metodo_newton_raphson(110.0)

r2_bis, iter_r2_bis = metodo_biseccion(2160.0, 2340.0)
r2_nw,  iter_r2_nw  = metodo_newton_raphson(2250.0)

# --- CÁLCULO EXPLÍCITO DE SENSIBILIDAD (VIÑETA 3) ---
# Evaluamos la derivada analítica exactamente en el valor de la segunda raíz hallada (~2250 Hz)
_, dZ_df_r2 = funcion_objetivo(r2_nw)
sensibilidad_física = 1.0 / dZ_df_r2

# --- REPORTE DE IMPRESIÓN COMPILADO Y LIMPIO ---
print("======================================================================================")
print("DETECCIÓN DE RAÍCES, COMPARACIÓN DE MÉTODOS Y SENSIBILIDAD")
print("======================================================================================")
print("  1: LÍMITES DE LA BANDA DE OPERACIÓN SEGURA (|Z| <= 150 \u03a9)")
print(f"  -> Límite Inferior (Raíz f_1):   {r1_nw:.4f} Hz  (Calculada con Newton-Raphson)")
print(f"  -> Límite Superior (Raíz f_2):  {r2_nw:.4f} Hz  (Calculada con Newton-Raphson)")
print(f"  >> RANGO DE OPERACIÓN SEGURA:   [{r1_nw:.2f} Hz , {r2_nw:.2f} Hz]")
print("--------------------------------------------------------------------------------------")
print("  2: COMPARACIÓN DE EFICIENCIA (NÚMERO DE ITERACIONES)")
print("   Objetivo de Búsqueda      |  Iteraciones BISECCIÓN  |  Iteraciones NEWTON-RAPHSON")
print("--------------------------------------------------------------------------------------")
print(f"   Raíz f_1 (~114 Hz)        |           {iter_r1_bis:2d}            |             {iter_r1_nw:2d}")
print(f"   Raíz f_2 (~2244 Hz)       |           {iter_r2_bis:2d}            |             {iter_r2_nw:2d}")
print("  ------------------------------------------------------------------------------------")
print("  3: ANÁLISIS DE SENSIBILIDAD EN LA RAÍZ SUPERIOR (f \u2248 2244 Hz)")
print(f"  -> Pendiente espectral en f_2 (d|Z|/df):   {dZ_df_r2:.6f} \u03a9/Hz")
print(f"  -> SENSIBILIDAD INVERSA (df/d|Z|):          {sensibilidad_física:.4f} Hz/\u03a9")
print("======================================================================================")

# --- SISTEMA GRÁFICO SEGURO DE VALIDACIÓN ---
try:
    import matplotlib.pyplot as plt
    f_eje = np.linspace(frecuencias.min(), frecuencias.max(), 1000)
    z_curva = [evaluar_spline_y_derivada(frecuencias, impedancias, M_global, h_global, fx)[0] for fx in f_eje]
    
    plt.figure(figsize=(9, 4.8))
    plt.plot(f_eje, z_curva, color='#1f77b4', linewidth=2.0, label=r'Perfil de Impedancia $|Z|(f)$')
    plt.axhline(150, color='#d62728', linestyle='--', linewidth=1.2, label=r'Umbral Crítico $Z_{th} = 150\ \Omega$')
    
    # Marcar de forma explícita los dos cruces por cero en el gráfico
    plt.scatter([r1_nw, r2_nw], [150, 150], color='black', s=50, zorder=5)
    plt.axvline(r1_nw, color='gray', linestyle=':', alpha=0.6)
    plt.axvline(r2_nw, color='gray', linestyle=':', alpha=0.6)
    
    plt.title('Identificación de la Banda Operativa Segura mediante Cruce de Umbral', fontsize=11, fontweight='bold', pad=12)
    plt.xlabel('Frecuencia $f$ (Hz)', fontsize=10)
    plt.ylabel(r'Impedancia $|Z|$ ($\Omega$)', fontsize=10)
    plt.xlim(50, 2850)
    plt.ylim(128, 162)
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.legend(loc='lower right', fontsize=9)
    
    plt.savefig('grafico_banda_segura_D.png', dpi=300, bbox_inches='tight')
    plt.show()
except ModuleNotFoundError:
    pass
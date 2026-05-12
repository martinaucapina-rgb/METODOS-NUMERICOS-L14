import numpy as np
from scipy.interpolate import lagrange

# ======================================================================
# EJERCICIO 2: DIFERENCIACIÓN CON DATOS TABULADOS (x=2.25)
# ======================================================================

# 1. Datos extraídos de la tabla proporcionada
x_tabla = np.array([1.5, 1.9, 2.1, 2.4, 2.6, 3.1])
y_tabla = np.array([1.0628, 1.3961, 1.5432, 1.7349, 1.8423, 2.0397])
x_eval = 2.25

# 2. PROCESAMIENTO MATEMÁTICO
# Construimos el polinomio de Lagrange
P = lagrange(x_tabla, y_tabla)

# Obtenemos las funciones derivadas analíticas
P_derivada_1 = np.polyder(P, 1)
P_derivada_2 = np.polyder(P, 2)

# Evaluamos en el punto solicitado
res_f1 = P_derivada_1(x_eval)
res_f2 = P_derivada_2(x_eval)

# --- ENCABEZADO Y DATOS INICIALES ---
print("="*100)
print("       EJERCICIO 2: ANÁLISIS DE DIFERENCIACIÓN CON DATOS TABULADOS")
print("="*100)
print(f" Punto de interés para f'(x) y f''(x): x = {x_eval}")

print("\nDATOS INICIALES (TABLA):")
print("-" * 45)
print(f"{'Punto (i)':<12} | {'x_i':<10} | {'f(x_i)':<15}")
print("-" * 45)
for i in range(len(x_tabla)):
    print(f"Punto {i+1:<7} | {x_tabla[i]:<10} | {y_tabla[i]:<15.4f}")
print("-" * 100)

# --- SECCIÓN DE METODOLOGÍA ---
print("\n>>> METODOLOGÍA: INTERPOLACIÓN Y DERIVACIÓN ANALÍTICA")
print(" Fórmulas base:")
print("  1. Polinomio de Lagrange: P(x) = Σ [ y_i * L_i(x) ]")
print("  2. Aproximación f'(x):    d/dx [P(x)]")
print("  3. Aproximación f''(x):   d²/dx² [P(x)]")
print("-" * 100)

# --- SECCIÓN DE PROCESO (PASO A PASO) ---
print("\nPASO 1: CONSTRUCCIÓN DE LA FUNCIÓN P(x)")
print(" Utilizando los 6 puntos, se obtiene el siguiente polinomio de grado 5:")
print(f"\n P(x) =\n{np.poly1d(P)}")
print("-" * 100)

print("\nPASO 2: DERIVACIÓN TÉRMINO A TÉRMINO")
print(f" Se calculan las derivadas de P(x) para evaluar en x = {x_eval}:")
print(f"\n f'(x) ≈ P'(x) =\n{np.poly1d(P_derivada_1)}")
print(f"\n f''(x) ≈ P''(x) =\n{np.poly1d(P_derivada_2)}")
print("-" * 100)

# --- RESULTADOS FINALES ---
print(f"\nPASO 3: RESULTADOS FINALES EN x = {x_eval}")
print("-" * 60)

print(f" PRIMERA DERIVADA:")
print(f"  f'({x_eval}) = {res_f1:.6f}")

print("-" * 30)

print(f" SEGUNDA DERIVADA:")
print(f"  f''({x_eval}) = {res_f2:.8f}")

print("\n" + "="*100)
print("            FIN DEL EJERCICIO 2")
print("="*100)

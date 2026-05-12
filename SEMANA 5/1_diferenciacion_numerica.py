import numpy as np

def f(x):
    return np.exp(x)

# Configuración del análisis
x_val = 2
v_real = np.exp(2)
h_valores = [0.5, 0.1, 0.05, 0.01, 0.005, 0.001]

def calcular_error(aprox, real):
    return abs((real - aprox) / real) * 100

# ======================================================================
# EJERCICIO 1: DIFERENCIACIÓN NUMÉRICA PARA f(x) = e^x
# ======================================================================
print("="*100)
print("       EJERCICIO 1: ANÁLISIS DE DIFERENCIACIÓN NUMÉRICA")
print("="*100)
print(f" Función: f(x) = e^x")
print(f" Punto de evaluación: x = {x_val}")
print(f" Valor real exacto (e^2): {v_real:.10f}")
print("="*100)

# --- SECCIÓN 1: PRIMERA DERIVADA ---
print(f"\n>>> 1.1 PRIMERA DERIVADA f'({x_val})")
print(" Fórmulas utilizadas:")
print("  - Regresiva:  [f(x) - f(x-h)] / h")
print("  - Progresiva: [f(x+h) - f(x)] / h")
print("  - Central:    [f(x+h) - f(x-h)] / 2h")
print("-" * 165)
print(f"{'h':<8} | {'Regresiva':<15} | {'Err% Reg':<12} | {'Progresiva':<15} | {'Err% Prog':<12} | {'Central':<15} | {'Err% Cent':<12}")
print("-" * 165)

for h in h_valores:
    f1_reg = (f(x_val) - f(x_val - h)) / h
    f1_pro = (f(x_val + h) - f(x_val)) / h
    f1_cen = (f(x_val + h) - f(x_val - h)) / (2 * h)
    
    e_reg = f"{calcular_error(f1_reg, v_real):.6f}%"
    e_pro = f"{calcular_error(f1_pro, v_real):.6f}%"
    e_cen = f"{calcular_error(f1_cen, v_real):.6f}%"
    
    # Impresión con 6 decimales para la primera derivada
    print(f"{h:<8} | {f1_reg:<15.6f} | {e_reg:<12} | {f1_pro:<15.6f} | {e_pro:<12} | {f1_cen:<15.6f} | {e_cen:<12}")

# --- SECCIÓN 2: SEGUNDA DERIVADA ---
print(f"\n>>> 1.2 SEGUNDA DERIVADA f''({x_val})")
print(" Fórmulas utilizadas:")
print("  - Regresiva:  [f(x) - 2f(x-h) + f(x-2h)] / h^2")
print("  - Progresiva: [f(x+2h) - 2f(x+h) + f(x)] / h^2")
print("  - Central:    [f(x+h) - 2f(x) + f(x-h)] / h^2")
print("-" * 185)
print(f"{'h':<8} | {'Regresiva':<18} | {'Err% Reg':<15} | {'Progresiva':<18} | {'Err% Prog':<15} | {'Central':<18} | {'Err% Cent':<15}")
print("-" * 185)

for h in h_valores:
    f2_reg = (f(x_val) - 2*f(x_val - h) + f(x_val - 2*h)) / (h**2)
    f2_pro = (f(x_val + 2*h) - 2*f(x_val + h) + f(x_val)) / (h**2)
    f2_cen = (f(x_val + h) - 2*f(x_val) + f(x_val - h)) / (h**2)
    
    e_reg2 = f"{calcular_error(f2_reg, v_real):.8f}%"
    e_pro2 = f"{calcular_error(f2_pro, v_real):.8f}%"
    e_cen2 = f"{calcular_error(f2_cen, v_real):.8f}%"
    
    # Impresión con 8 decimales para la segunda derivada
    print(f"{h:<8} | {f2_reg:<18.8f} | {e_reg2:<15} | {f2_pro:<18.8f} | {e_pro2:<15} | {f2_cen:<18.8f} | {e_cen2:<15}")

print("\n" + "="*100)
print("            FIN DEL EJERCICIO 1")
print("="*100)

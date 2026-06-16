import numpy as np
import matplotlib.pyplot as plt

# =================================================================
# 1. DEFINICIÓN DEL OSCILADOR DE VAN DER POL (SISTEMA ODE)
# =================================================================
def van_der_pol(t, estado):
    x, v = estado
    mu = 1.5
    dxdt = v
    dvdt = mu * (1 - x**2) * v - x
    return np.array([dxdt, dvdt])

# =================================================================
# 2. PARÁMETROS DE SIMULACIÓN
# =================================================================
t_0 = 0.0
t_f = 20.0
h = 0.05
Num_pasos = int((t_f - t_0) / h)

estado_inicial = np.array([2.0, 0.0]) # x(0) = 2.0, v(0) = 0.0

historial_t = np.zeros(Num_pasos + 1)
historial_estado = np.zeros((Num_pasos + 1, 2))

historial_t[0] = t_0
historial_estado[0] = estado_inicial

# =================================================================
# 3. ALGORITMO RUNGE-KUTTA DE 4TO ORDEN (RK4)
# =================================================================
estado_actual = estado_inicial
t_actual = t_0

for i in range(Num_pasos):
    k1 = van_der_pol(t_actual, estado_actual)
    k2 = van_der_pol(t_actual + h/2, estado_actual + (h/2) * k1)
    k3 = van_der_pol(t_actual + h/2, estado_actual + (h/2) * k2)
    k4 = van_der_pol(t_actual + h, estado_actual + h * k3)
    
    estado_actual = estado_actual + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
    t_actual += h
    
    historial_t[i+1] = t_actual
    historial_estado[i+1] = estado_actual

x_vals = historial_estado[:, 0]
v_vals = historial_estado[:, 1]

# =================================================================
# 4. TABLA DE DATOS EN CONSOLA (ENTREGABLE OBLIGATORIO)
# =================================================================
print("="*68)
print("          TABLA DE RESULTADOS SIMULACIÓN RK4 - VAN DER POL")
print("="*68)
print(f" {'Iteración':^9} | {'Tiempo (s)':^12} | {'Posición x(t)':^18} | {'Velocidad v(t)':^18}")
print("-"*68)

for i in range(11):
    print(f" {i:^9} | {historial_t[i]:^12.2f} | {x_vals[i]:^18.6f} | {v_vals[i]:^18.6f}")

print(f" {'...':^9} | {'...':^12} | {'...':^18} | {'...':^18}")
print(f" {Num_pasos:^9} | {historial_t[-1]:^12.2f} | {x_vals[-1]:^18.6f} | {v_vals[-1]:^18.6f}")
print("="*68)

# =================================================================
# 5. GENERACIÓN Y GUARDADO DE GRÁFICAS (ENTREGABLES OBLIGATORIOS)
# =================================================================
# Gráfica 1: Evolución temporal
plt.figure(figsize=(10, 4))
plt.plot(historial_t, x_vals, label="Desplazamiento x(t)", color="royalblue", lw=2)
plt.plot(historial_t, v_vals, label="Velocidad v(t)", color="crimson", linestyle="--", lw=1.5)
plt.title("Gráfica 1: Evolución Temporal del Desplazamiento y Velocidad (RK4)")
plt.xlabel("Tiempo t (segundos)")
plt.ylabel("Amplitud")
plt.grid(True, linestyle=":", alpha=0.6)
plt.xlim(t_0, t_f)
plt.legend(loc="upper right")

# Guardar Gráfica 1 antes de mostrarla (dpi=300 asegura alta calidad para el informe)
plt.savefig("grafica_1_evolucion_temporal.png", dpi=300, bbox_inches='tight')
plt.show()

# Gráfica 2: Espacio de fases
plt.figure(figsize=(5, 5))
plt.plot(x_vals, v_vals, color="darkviolet", lw=2)
plt.title("Gráfica 2: Espacio de Fases (Ciclo Límite)")
plt.xlabel("Desplazamiento x(t)")
plt.ylabel("Velocidad v(t)")
plt.grid(True, linestyle=":", alpha=0.6)
plt.tight_layout()

# Guardar Gráfica 2 antes de mostrarla
plt.savefig("grafica_2_espacio_fases.png", dpi=300, bbox_inches='tight')
plt.show()
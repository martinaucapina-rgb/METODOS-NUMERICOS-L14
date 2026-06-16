import numpy as np

# 1. DEFINICIÓN DEL OSCILADOR DE VAN DER POL
def van_der_pol(t, estado):
    x, v = estado
    mu = 1.5
    dxdt = v
    dvdt = mu * (1 - x**2) * v - x
    return np.array([dxdt, dvdt])

# 2. PARÁMETROS DE SIMULACIÓN (PASO h = 0.5)
t_0 = 0.0
t_f = 20.0
h = 0.5  # <--- Incremento drástico del paso
Num_pasos = int((t_f - t_0) / h)

estado_inicial = np.array([2.0, 0.0])

historial_t = np.zeros(Num_pasos + 1)
historial_estado = np.zeros((Num_pasos + 1, 2))

historial_t[0] = t_0
historial_estado[0] = estado_inicial

# 3. ALGORITMO RK4
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

# 4. IMPRESIÓN DE LA TABLA EN CONSOLA
print("="*68)
print("          TABLA DE RESULTADOS SIMULACIÓN RK4 (h = 0.5 s)")
print("="*68)
print(f" {'Iteración':^9} | {'Tiempo (s)':^12} | {'Posición x(t)':^18} | {'Velocidad v(t)':^18}")
print("-"*68)

# Mostramos las primeras 11 iteraciones para ver la divergencia inicial
for i in range(min(11, Num_pasos + 1)):
    print(f" {i:^9} | {historial_t[i]:^12.2f} | {x_vals[i]:^18.6e} | {v_vals[i]:^18.6e}")
print("="*68)
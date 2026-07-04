import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# ALGORITMO: Solución de EDOs (Método de Runge-Kutta de 4to Orden - RK4)
# =============================================================================
def ejecutar_simulacion_rk4():
    """
    Simula dinámicamente el transitorio de voltaje empleando parámetros sintonizados
    (R*C = 895.3) para lograr un acoplamiento óptimo con la Tabla 1.
    """
    print("=====================================================================")
    print("         SOLUCIÓN DE EDOS: RUNGE-KUTTA DE 4TO ORDEN (RK4)           ")
    print("=====================================================================")
    
    # Parámetros sintonizados para calzar exactamente con la curva experimental
    V_asintota = 0.0   # Nivel de referencia de descarga completa (V)
    R = 13.90          # Resistencia interna sintonizada (Ohms)
    C = 64.41          # Capacitancia equivalente sintonizada (F)
    tau = R * C        # Constante de tiempo del circuito (tau = 895.3 s)
    
    # Configuración del espacio de integración temporal
    t_0 = 0.0          
    t_final = 60.0     
    h = 2.0            # Tamaño de paso constante (s)
    
    # Condición inicial extraída estrictamente de la Tabla 1
    V_0 = 48.00        
    
    # EDO despejada: dV/dt = f(t, V)
    def f(t, V):
        return (V_asintota - V) / tau
    
    t_actual = t_0
    V_actual = V_0
    
    historial_tiempo = [t_actual]
    historial_voltaje = [V_actual]
    
    print("-> Parámetros de la simulación Sintonizada Dinámicamente:")
    print(f"   EDO: dV/dt = ({V_asintota} - V) / {tau:.1f}")
    print(f"   Intervalo: [{t_0}, {t_final}] s con h = {h} s")
    print(f"   Condición inicial: V({t_0}) = {V_0} V\n")
    
    print("---------------------------------------------------------------------")
    print("  PASO (i)  |   TIEMPO (s)  |    k1   |    k2   |    k3   | VOLTAJE (V)")
    print("---------------------------------------------------------------------")
    print(f"    00      |     {t_actual:4.1f}      |    --   |    --   |    --   |   {V_actual:6.4f}")
    
    paso = 0
    while t_actual < t_final:
        paso += 1
        
        # Cómputo de las 4 pendientes del subintervalo (RK4)
        k1 = f(t_actual, V_actual)
        k2 = f(t_actual + h/2.0, V_actual + (h/2.0)*k1)
        k3 = f(t_actual + h/2.0, V_actual + (h/2.0)*k2)
        k4 = f(t_actual + h, V_actual + h*k3)
        
        # Ecuación recursiva de actualización de Runge-Kutta 4
        V_siguiente = V_actual + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        t_siguiente = t_actual + h
        
        f_k1 = f"{k1:.4f}"
        f_k2 = f"{k2:.4f}"
        f_k3 = f"{k3:.4f}"
        f_volt = f"{V_siguiente:.2f}" if V_siguiente == round(V_siguiente, 2) else f"{V_siguiente:.4f}"
        
        if paso <= 3 or t_siguiente >= t_final - 2.0:
            print(f"    {paso:02d}      |     {t_siguiente:4.1f}      | {f_k1} | {f_k2} | {f_k3} |   {f_volt}")
        elif paso == 4:
            print("    ...     |      ...      |   ...   |   ...   |   ...   |     ...")
            
        t_actual = t_siguiente
        V_actual = V_siguiente
        
        historial_tiempo.append(t_actual)
        historial_voltaje.append(V_actual)
        
    print("---------------------------------------------------------------------")
    print(f"RESULTADO: Simulación completada con éxito ({paso} pasos).")
    print(f"Voltaje teórico sintonizado final: V({t_actual:.1f}) = {V_actual:.4f} V.")
    print("=====================================================================\n")
    
    return np.array(historial_tiempo), np.array(historial_voltaje)


# =============================================================================
# BLOQUE DE VALIDACIÓN GRÁFICA OPTIMIZADA
# =============================================================================
if __name__ == "__main__":
    # 1. Obtener curva continua de la simulación numérica RK4 sintonizada
    t_sim, v_sim = ejecutar_simulacion_rk4()
    
    # 2. Transcripción exacta de las lecturas experimentales reales (Tabla 1)
    tiempo_real = np.arange(0.0, 65.0, 5.0)
    voltaje_real = np.array([48.00, 47.85, 47.60, 47.32, 47.05, 46.80, 46.52, 
                             46.20, 45.90, 45.65, 45.30, 45.02, 44.75])
    
    # 3. Construcción del gráfico comparativo de alta fidelidad
    plt.figure(figsize=(10, 6))
    
    # Curva continua para el modelo RK4 corregido
    plt.plot(t_sim, v_sim, color='#1f77b4', linestyle='-', linewidth=2.5, 
             label='Modelo Teórico Transitorio Sintonizado (RK4)')
    
    # Marcadores para las mediciones empíricas reales de la Tabla 1
    plt.scatter(tiempo_real, voltaje_real, color='#d62728', marker='o', s=60, zorder=5,
                label='Voltaje de Línea Real Vdc (Tabla 1)')
    
    # Personalización estética
    plt.title('Validación de la Dinámica de Voltaje: Circuito Sintonizado vs Datos Reales', 
              fontsize=12, fontweight='bold', pad=15)
    plt.xlabel('Tiempo Transcurrido (s)', fontsize=11, fontweight='bold')
    plt.ylabel('Tensión del Sistema Vdc (V)', fontsize=11, fontweight='bold')
    
    plt.xlim(-2, 62)
    plt.ylim(44.0, 48.5)
    
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc='upper right', fontsize=10.5, frameon=True, shadow=True)
    
    plt.tight_layout()
    plt.show()
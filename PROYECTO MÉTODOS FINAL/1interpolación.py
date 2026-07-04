import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

# =============================================================================
# BLOQUE GLOBAL: Base de Datos Experimentales (Muestreo a 60 segundos)
# =============================================================================
# Vectores numéricos que indexan de forma cronológica los datos de la Tabla 1
tiempo = np.arange(0.0, 65.0, 5.0)  # Variable independiente (s)
v_dc   = np.array([48.00, 47.85, 47.60, 47.32, 47.05, 46.80, 46.52, 
                   46.20, 45.90, 45.65, 45.30, 45.02, 44.75])  # Tensión (V)
i_fase = np.array([0.0, 12.5, 18.2, 22.0, 24.5, 26.1, 27.0, 
                   27.5, 27.8, 27.9, 28.0, 27.6, 26.8])        # Corriente (A)
t_est  = np.array([25.0, 25.2, 25.6, 26.1, 26.8, 27.5, 28.3, 
                   29.2, 30.1, 31.0, 32.0, 33.1, 34.2])        # Temperatura (°C)


# =============================================================================
# ALGORITMO: Interpolación de Tensión mediante Splines Cúbicos Naturales
# =============================================================================
def ejecutar_interpolacion_trazadores():
    """
    Desarrolla la interpolación por trazadores cúbicos seleccionando el subconjunto
    de datos que encierran el punto de interés, evaluando las condiciones de contorno.
    """
    print("=====================================================================")
    print("      PROCESAMIENTO DIGITAL DE INTERPOLACIÓN: SPLINES CÚBICOS        ")
    print("=====================================================================")
    
    # Segmentación exacta del intervalo de análisis alrededor de t = 27.5 s
    # Se extraen los 4 nodos consecutivos requeridos para el planteamiento matricial
    t_nodos = np.array([tiempo[3], tiempo[4], tiempo[5], tiempo[6]])  # [15.0, 20.0, 25.0, 30.0]
    v_nodos = np.array([v_dc[3], v_dc[4], v_dc[5], v_dc[6]])      # [47.32, 47.05, 46.80, 46.52]
    
    print(f"-> Nodo Inicial t_0 = {t_nodos[0]} s | Tensión V(t_0) = {v_nodos[0]} V")
    print(f"-> Nodo Interno t_1 = {t_nodos[1]} s | Tensión V(t_1) = {v_nodos[1]} V")
    print(f"-> Nodo Interno t_2 = {t_nodos[2]} s | Tensión V(t_2) = {v_nodos[2]} V")
    print(f"-> Nodo Terminal t_3 = {t_nodos[3]} s | Tensión V(t_3) = {v_nodos[3]} V\n")
    
    # Configuración matemática del Spline Cúbico Natural
    # bc_type='natural' impone analíticamente que la curvatura en las fronteras sea nula
    trazador_sistema = CubicSpline(t_nodos, v_nodos, bc_type='natural')
    
    # Extracción de las segundas derivadas para corroborar el desarrollo algebraico manual
    # La función derivative(2) evalúa la aceleración del cambio de voltaje en cada nodo
    d2_t0 = trazador_sistema.derivative(2)(t_nodos[0])
    d2_t1 = trazador_sistema.derivative(2)(t_nodos[1])
    d2_t2 = trazador_sistema.derivative(2)(t_nodos[2])
    d2_t3 = trazador_sistema.derivative(2)(t_nodos[3])
    
    print("---------------------------------------------------------------------")
    print("VALIDACIÓN DE SEGUNDAS DERIVADAS (MOMENTOS DE MATRIZ DE ACOPLAMIENTO)")
    print("---------------------------------------------------------------------")
    print(f"V''(t_0 = {t_nodos[0]}s) = {d2_t0:.5f} (Frontera Natural)")
    print(f"V''(t_1 = {t_nodos[1]}s) = {d2_t1:.5f} (Coincidencia con cálculo manual)")
    print(f"V''(t_2 = {t_nodos[2]}s) = {d2_t2:.5f} (Coincidencia con cálculo manual)")
    print(f"V''(t_3 = {t_nodos[3]}s) = {d2_t3:.5f} (Frontera Natural)\n")
    
    # Evaluación numérica en el punto intermedio no tabulado
    t_objetivo = 27.5
    v_interp = trazador_sistema(t_objetivo)
    
    print("---------------------------------------------------------------------")
    print("RESULTADO DE LA ESTIMACIÓN INTERMEDIA")
    print("---------------------------------------------------------------------")
    print(f"Para un instante de tiempo t = {t_objetivo} s")
    print(f"La tensión interpolada en el bus DC es: V_dc({t_objetivo}) = {v_interp:.4f} V")
    print("=====================================================================\n")

# Ejecución del procedimiento analítico
ejecutar_interpolacion_trazadores()
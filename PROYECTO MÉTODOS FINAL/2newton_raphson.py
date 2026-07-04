import numpy as np

# =============================================================================
# ALGORITMO INTEGRADO: MÍNIMOS CUADRADOS (CRAMER) Y BÚSQUEDA DE RAÍCES
# =============================================================================
def ejecutar_sistema_raices_completo():
    """
    Bloque automatizado que procesa el sistema de ecuaciones normales de tu reporte
    mediante la Regla de Cramer para deducir el polinomio térmico y validar las 
    iteraciones de Newton-Raphson de forma dinámica.
    """
    print("=====================================================================")
    print("      SISTEMA DE RAÍCES: AJUSTE POLINOMIAL Y NEWTON-RAPHSON          ")
    print("=====================================================================")
    
    # -------------------------------------------------------------------------
    # PASO 1: CONSTRUCCIÓN DEL SISTEMA DE ECUACIONES NORMALES (TABLA 1)
    # -------------------------------------------------------------------------
    # Se inicializan los coeficientes de la matriz con las sumatorias del informe
    A = np.array([
        [13.0,   390.0,    16250.0],
        [390.0,  16250.0,  760500.0],
        [16250.0, 760500.0, 37781250.0]
    ])
    B = np.array([376.5, 11721.5, 496667.5])
    
    # -------------------------------------------------------------------------
    # PASO 2: RESOLUCIÓN DEL AJUSTE MEDIANTE LA REGLA DE CRAMER
    # -------------------------------------------------------------------------
    # Cómputo del determinante de la matriz de coeficientes principal
    det_principal = np.linalg.det(A)
    
    # Determinante asociado al término independiente (a0)
    A0 = A.copy(); A0[:, 0] = B
    det_A0 = np.linalg.det(A0)
    
    # Determinante asociado al término lineal (a1)
    A1 = A.copy(); A1[:, 1] = B
    det_A1 = np.linalg.det(A1)
    
    # Determinante asociado al término cuadrático (a2)
    A2 = A.copy(); A2[:, 2] = B
    det_A2 = np.linalg.det(A2)
    
    # Cálculo dinámico de los coeficientes del polinomio continuo
    # Se formatean para acoplarse con la resolución manual de determinantes
    a0 = round(det_A0 / det_principal, 4)
    a1 = round(det_A1 / det_principal, 4)
    a2 = round(det_A2 / det_principal, 4)
    
    print("-> Coeficientes calculados mediante la Regla de Cramer:")
    print(f"   a0 (T_inicial) = {a0} °C")
    print(f"   a1 (Lin_Joule) = {a1}")
    print(f"   a2 (No_lineal) = {a2}")
    print(f"   Ecuación deducida: Te(t) = {a0} + {a1}t + {a2}t^2\n")
    
    # -------------------------------------------------------------------------
    # PASO 3: CONFIGURACIÓN DEL ALGORITMO DE NEWTON-RAPHSON
    # -------------------------------------------------------------------------
    T_critico = 32.0  # Umbral de seguridad del estator (°C)
    
    # Definición de la Función Objetivo: f(t) = Te(t) - T_critico = 0
    def f(t):
        return a2 * (t**2) + a1 * t + (a0 - T_critico)
    
    # Primera derivada analítica del modelo: f'(t) = 2*a2*t + a1
    def df(t):
        return 2 * a2 * t + a1
    
    # Parámetros del espacio de búsqueda numérico
    t_actual = 45.0          # Aproximación inicial t0 (s)
    tolerancia = 1e-3        # Umbral de parada (0.001%)
    max_iteraciones = 10     
    
    print("---------------------------------------------------------------------")
    print(" ITERACIÓN | TIEMPO ESTIMADO (s) |    f(t)   |   f'(t)   |  ERROR (%)")
    print("---------------------------------------------------------------------")
    
    paso = 0
    error_porcentual = 100.0
    
    while error_porcentual > tolerancia and paso < max_iteraciones:
        paso += 1
        
        # Evaluación del estado actual de pendientes
        f_val = f(t_actual)
        df_val = df(t_actual)
        
        # Proyección geométrica de la raíz (Newton-Raphson)
        t_siguiente = t_actual - (f_val / df_val)
        
        # Cálculo del error relativo aproximado porcentual
        error_porcentual = abs((t_siguiente - t_actual) / t_siguiente) * 100
        
        # Despliegue de la traza de ejecución del algoritmo
        print(f"    {paso:02d}     |      {t_siguiente:9.4f}      | {f_val:9.4f} | {df_val:8.4f}  | {error_porcentual:8.4f}%")
        
        # Transición de estado para el siguiente ciclo
        t_actual = t_siguiente
        
    print("---------------------------------------------------------------------")
    print(f"RESULTADO: Convergencia lograda con éxito en la iteración {paso}.")
    print(f"El umbral de alerta preventiva se alcanza a los t = {t_actual:.4f} s.")
    print("=====================================================================\n")

if __name__ == "__main__":
    ejecutar_sistema_raices_completo()
import numpy as np

# =============================================================================
# ALGORITMO: Derivación Numérica (Diferencias Finitas de Alta Precisión)
# =============================================================================
def ejecutar_derivacion_numerica(t_vector, i_vector, t_objetivo=20.0):
    """
    Calcula la primera derivada de la corriente mediante esquemas Forward, 
    Backward y Central de orden O(h^2) para un punto de muestreo dado.
    """
    print("=====================================================================")
    print("      DERIVACIÓN NUMÉRICA: DIFERENCIAS FINITAS DE ORDEN O(h^2)       ")
    print("=====================================================================")
    
    # Identificación automática del tamaño de paso 'h' (Muestreo homogéneo)
    h = t_vector[1] - t_vector[0]
    
    # Localizar el índice del nodo de interés en el vector de tiempo
    idx = np.where(t_vector == t_objetivo)[0][0]
    
    print(f"-> Analizando nodo en t = {t_objetivo} s (Índice de registro: {idx})")
    print(f"   Tamaño de paso constante (h) = {h} s\n")
    
    # 1. Esquema Hacia Adelante (Forward) de Orden O(h^2)
    derivada_forward = (-i_vector[idx + 2] + 4 * i_vector[idx + 1] - 3 * i_vector[idx]) / (2 * h)
    
    # 2. Esquema Hacia Atrás (Backward) de Orden O(h^2)
    derivada_backward = (3 * i_vector[idx] - 4 * i_vector[idx - 1] + i_vector[idx - 2]) / (2 * h)
    
    # 3. Esquema Central (Central) de Orden O(h^2)
    derivada_central = (i_vector[idx + 1] - i_vector[idx - 1]) / (2 * h)
    
    # Función auxiliar para formatear los resultados sin ceros innecesarios
    def formatear_resultado(valor):
        return f"{valor:.1f}" if valor == round(valor, 1) else f"{valor:.4f}"
    
    print("---------------------------------------------------------------------")
    print(" ESQUEMA APLICADO (Orden O(h^2))            |   RESULTADO (A/s)     ")
    print("---------------------------------------------------------------------")
    print(f" Diferencias Finitas Hacia Adelante (Forward) |        {formatear_resultado(derivada_forward)}")
    print(f" Diferencias Finitas Hacia Atrás (Backward)   |        {formatear_resultado(derivada_backward)}")
    print(f" Diferencias Finitas Centrales (Central)     |        {formatear_resultado(derivada_central)}")
    print("---------------------------------------------------------------------")
    print(f"RESULTADO FINAL: Tasa central óptima dIf/dt = {formatear_resultado(derivada_central)} A/s.")
    print("=====================================================================\n")

# =============================================================================
# BLOQUE DE PRUEBA: Integración con la Base de Datos de la Tabla 1
# =============================================================================
if __name__ == "__main__":
    # Vectores base del experimento de la moto eléctrica
    tiempo_eje = np.arange(0.0, 65.0, 5.0)
    i_fase_eje = np.array([0.0, 12.5, 18.2, 22.0, 24.5, 26.1, 27.0, 
                           27.5, 27.8, 27.9, 28.0, 27.6, 26.8])
    
    # Ejecución del módulo de cálculo para t = 20.0 segundos
    ejecutar_derivacion_numerica(tiempo_eje, i_fase_eje, t_objetivo=20.0)
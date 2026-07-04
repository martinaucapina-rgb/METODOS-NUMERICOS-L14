import numpy as np

# =============================================================================
# ALGORITMO: Integración Numérica (Regla de Simpson 1/3 Compuesta o Múltiple)
# =============================================================================
def ejecutar_integracion_simpson(t_vector, i_vector):
    """
    Calcula el área bajo la curva de la corriente (Carga eléctrica total Q)
    utilizando la Regla de Simpson 1/3 con n subintervalos (n debe ser par).
    """
    print("=====================================================================")
    print("      INTEGRACIÓN NUMÉRICA: REGLA DE SIMPSON 1/3 MÚLTIPLE            ")
    print("=====================================================================")
    
    n = len(t_vector) - 1  # Número de subintervalos
    h = t_vector[1] - t_vector[0]  # Tamaño de paso constante
    
    print(f"-> Parámetros del espacio de integración:")
    print(f"   Límites de integración: [a, b] = [{t_vector[0]}, {t_vector[-1]}] s")
    print(f"   Número de subintervalos (n) = {n} ({'Válido: Par' if n % 2 == 0 else 'Inválido: Impar'})")
    print(f"   Ancho de cada intervalo (h) = {h} s\n")
    
    if n % 2 != 0:
        print("[ERROR] El método de Simpson 1/3 compuesto requiere un número par de subintervalos.")
        return None

    # Extracción de componentes según el esquema de Simpson 1/3
    extremos = i_vector[0] + i_vector[-1]
    
    # Índices impares: 1, 3, 5, ..., n-1
    sum_impares = np.sum(i_vector[1:-1:2])
    
    # Índices pares internos: 2, 4, 6, ..., n-2
    sum_pares = np.sum(i_vector[2:-1:2])
    
    # Aplicación formal de la fórmula compuesta
    integral_q = (h / 3) * (extremos + 4 * sum_impares + 2 * sum_pares)
    
    # Función auxiliar para formatear los resultados sin ceros innecesarios
    def formatear_resultado(valor):
        return f"{valor:.1f}" if valor == round(valor, 1) else f"{valor:.4f}"
    
    print("---------------------------------------------------------------------")
    print(" COMPONENTE DEL ALGORITMO                   |   VALOR CONSOLIDADO    ")
    print("---------------------------------------------------------------------")
    print(f" Suma de Nodos Extremos [I(0) + I(60)]       |        {formatear_resultado(extremos)}")
    print(f" Suma de Nodos Impares (Ponderador 4)       |        {formatear_resultado(sum_impares)}")
    print(f" Suma de Nodos Pares Internos (Ponderador 2)|        {formatear_resultado(sum_pares)}")
    print("---------------------------------------------------------------------")
    print(f"RESULTADO FINAL: Carga total calculada Q = {formatear_resultado(integral_q)} C.")
    print("=====================================================================\n")

# =============================================================================
# BLOQUE DE PRUEBA: Integración con la Base de Datos de la Tabla 1
# =============================================================================
if __name__ == "__main__":
    # Vectores base del experimento de la moto eléctrica
    tiempo_eje = np.arange(0.0, 65.0, 5.0)
    i_fase_eje = np.array([0.0, 12.5, 18.2, 22.0, 24.5, 26.1, 27.0, 
                           27.5, 27.8, 27.9, 28.0, 27.6, 26.8])
    
    # Ejecución del módulo de integración
    ejecutar_integracion_simpson(tiempo_eje, i_fase_eje)
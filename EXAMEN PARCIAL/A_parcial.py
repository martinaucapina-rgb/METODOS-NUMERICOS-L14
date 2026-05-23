import numpy as np

# Definición de vectores en memoria para nodos independientes y dependientes
frecuencias = np.array([
    100, 120, 145, 170, 200, 235, 270, 310, 355, 405, 460, 520, 585, 655, 730,
    810, 895, 985, 1080, 1180, 1290, 1410, 1540, 1680, 1830, 1990, 2160, 2340, 2530, 2730
], dtype=float)

impedancias = np.array([
    152.3, 149.1, 146.8, 144.9, 142.0, 139.5, 137.9, 136.1, 134.8, 133.6, 132.7, 131.9, 131.4, 131.1, 130.9,
    131.0, 131.3, 131.9, 132.7, 133.8, 135.2, 136.9, 138.9, 141.1, 143.5, 146.1, 149.0, 152.2, 155.6, 159.2
], dtype=float)

# --- CÁLCULOS PRINCIPALES ---
# Para identificar el punto más bajo de los datos
idx_minimo = np.argmin(impedancias)
f_min_discreta = frecuencias[idx_minimo]
z_min_discreta = impedancias[idx_minimo]

print("--------------------------------------------------")
print("             ANÁLISIS EXPLORATORIO")
print("--------------------------------------------------")
print(f"Punto mínimo identificado en los datos tabulados:")
print(f"  - Índice de medición (i): {idx_minimo + 1}")
print(f"  - Frecuencia medida (f):  {f_min_discreta:.4f} Hz")
print(f"  - Impedancia medida (|Z|): {z_min_discreta:.4f} \u03a9")
print("--------------------------------------------------")

# --- GRAFICACIÓN ---
try:  # Para que no se rompa si se ejecuta en compiladores web sin librerías gráficas
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(7.5, 4.8))
    plt.plot(frecuencias, impedancias, color='#d62728', marker='o', 
             linestyle='-', linewidth=1.2, markersize=4.5, label='Mediciones Experimentales')
    
    plt.axvline(x=f_min_discreta, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    plt.axhline(y=z_min_discreta, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    
    plt.plot(f_min_discreta, z_min_discreta, color='#1f77b4', marker='x', 
             markersize=9, markeredgewidth=2, label=f'Mínimo Estimado \u2248 {f_min_discreta:.1f} Hz')
    
    plt.title('Perfil Espectroscópico: Magnitud de Impedancia $|Z|$ vs. Frecuencia $f$', 
              fontsize=11, fontweight='bold', pad=12)
    plt.xlabel('Frecuencia $f$ (Hz)', fontsize=9.5)
    plt.ylabel(r'Magnitud de Impedancia $|Z|$ ($\Omega$)', fontsize=9.5)
    plt.xlim(0, 2900)
    plt.ylim(125, 165)
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.legend(loc='upper center', fontsize=9, frameon=True)
    
    plt.savefig('grafico_parte_A.png', dpi=300, bbox_inches='tight')
    plt.show()

except ModuleNotFoundError:
    print("\n[INFO] Gráfica omitida automáticamente por restricciones del entorno en línea.")
    print("-----------------------------------------------------------------")
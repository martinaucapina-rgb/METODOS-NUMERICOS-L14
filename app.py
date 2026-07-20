import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

# ==========================================
# CONFIGURACIÓN DE LA PÁGINA
# ==========================================
st.set_page_config(
    page_title="Moto Eléctrica - Métodos Numéricos UNMSM",
    page_icon="⚡",
    layout="wide"
)

# TÍTULO PRINCIPAL
st.title("⚡ Dashboard de Métodos Numéricos: Tren de Potencia Síncrono")
st.subheader("Análisis Térmico, Electrónico y Energético de Motocicleta Eléctrica")
st.markdown("---")

# ==========================================
# BASE DE DATOS EXPERIMENTALES (TABLA 1)
# ==========================================
t_exp = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], dtype=float)
v_exp = np.array([48.00, 47.85, 47.60, 47.32, 47.05, 46.80, 46.52, 46.20, 45.85, 45.50, 45.12, 44.90, 44.75], dtype=float)
i_exp = np.array([0.0, 28.0, 27.5, 27.0, 26.8, 26.5, 26.2, 26.0, 25.8, 25.5, 25.3, 25.1, 25.0], dtype=float)
temp_exp = np.array([26.4, 26.7, 27.1, 27.5, 28.0, 28.6, 29.2, 29.9, 30.6, 31.3, 32.1, 32.9, 33.8], dtype=float)

# PANEL LATERAL
st.sidebar.header("⚙️ Panel de Control")
st.sidebar.subheader("Datos Medidos en Banco de Pruebas")

if st.sidebar.checkbox("Mostrar Tabla 1 (Datos Físicos)", value=True):
    df_tabla = pd.DataFrame({
        "Tiempo t (s)": t_exp,
        "Voltaje Vdc (V)": v_exp,
        "Corriente Im (A)": i_exp,
        "Temp. Estator Te (°C)": temp_exp
    })
    st.sidebar.dataframe(df_tabla, use_container_width=True)

# ESTRUCTURA DE 5 PESTAÑAS PARA CADA MÉTODO
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "1. Interpolación (Splines)",
    "2. Mínimos Cuadrados",
    "3. Newton-Raphson",
    "4. Derivación Numérica",
    "5. EDO (RK4) e Integración"
])

# ==========================================
# MÉTODO 1: INTERPOLACIÓN POR SPLINES CÚBICOS
# ==========================================
with tab1:
    st.header("1. Interpolación por Splines Cúbicos Naturales")
    st.write("Estimación no tabulada de la tensión instantánea en el bus de corriente continua ($V_{dc}$).")
    
    col1_1, col1_2 = st.columns([1, 2])
    
    with col1_1:
        st.subheader("Evaluación Puntual")
        t_eval = st.number_input("Seleccionar tiempo a evaluar t (s):", min_value=0.0, max_value=60.0, value=27.5, step=0.5)
        
        # Spline cúbico natural sobre los 4 nodos específicos del informe (t=15, 20, 25, 30)
        t_sub = np.array([15.0, 20.0, 25.0, 30.0])
        v_sub = np.array([47.32, 47.05, 46.80, 46.52])
        cs_sub = CubicSpline(t_sub, v_sub, bc_type='natural')
        v_calc_sub = float(cs_sub(t_eval))
        
        st.metric(label=f"Tensión estimada en t = {t_eval} s (Tramo Local)", value=f"{v_calc_sub:.4f} V")
        st.info(f"📌 Valor exacto del informe analítico manual para t = 27.5 s: **46.6635 V**")
        
        st.markdown("""
        **Segundas Derivadas calculadas:**
        * $V''(t_0 = 15s) = 0.00000$ (Frontera)
        * $V''(t_1 = 20s) = 0.00176$
        * $V''(t_2 = 25s) = -0.00224$
        * $V''(t_3 = 30s) = 0.00000$ (Frontera)
        """)
        
    with col1_2:
        st.subheader("Gráfica del Trazador Cúbico")
        t_dense = np.linspace(15, 30, 100)
        v_dense = cs_sub(t_dense)
        
        fig1, ax1 = plt.subplots(figsize=(8, 4))
        ax1.plot(t_sub, v_sub, 'ro', markersize=8, label='Nodos Experimentales')
        ax1.plot(t_dense, v_dense, 'b-', label='Spline Cúbico Natural')
        ax1.plot(t_eval, v_calc_sub, 'go', markersize=10, label=f'Punto Interpolado ({t_eval}s, {v_calc_sub:.4f}V)')
        ax1.set_xlabel("Tiempo t (s)")
        ax1.set_ylabel("Voltaje Vdc (V)")
        ax1.set_title("Interpolación en el Intervalo [15s, 30s]")
        ax1.grid(True, linestyle='--')
        ax1.legend()
        st.pyplot(fig1)

# ==========================================
# MÉTODO 2: MÍNIMOS CUADRADOS
# ==========================================
with tab2:
    st.header("2. Ajuste de Curvas por Mínimos Cuadrados Cuadráticos")
    st.write("Modelado continuo del calentamiento del estator mediante la resolución de la Matriz de Cramer.")
    
    a0, a1, a2 = 26.378, 0.0688, 0.0004
    
    col2_1, col2_2 = st.columns([1, 2])
    
    with col2_1:
        st.subheader("Ecuación Ajustada")
        st.latex(rf"T_e(t) = {a0:.4f} + {a1:.4f}t + {a2:.4f}t^2")
        
        st.markdown("""
        **Estructura del Sistema Lineal $M \\cdot A = B$:**
        * $\sum 1 = 13$
        * $\sum t = 390$
        * $\sum t^2 = 16250$
        * $\sum t^3 = 760500$
        * $\sum t^4 = 37927500$
        """)
        
        t_test = st.slider("Probar tiempo t (s):", 0.0, 60.0, 30.0, 1.0)
        temp_calc = a0 + a1 * t_test + a2 * (t_test**2)
        st.success(f"Temperatura estimada a los {t_test} s: **{temp_calc:.4f} °C**")

    with col2_2:
        st.subheader("Ajuste del Modelo vs. Mediciones Reales")
        t_model = np.linspace(0, 60, 200)
        temp_model = a0 + a1 * t_model + a2 * (t_model**2)
        
        fig2, ax2 = plt.subplots(figsize=(8, 4))
        ax2.plot(t_exp, temp_exp, 'rs', label='Lecturas de Sensores')
        ax2.plot(t_model, temp_model, 'g-', linewidth=2, label='Polinomio Ajustado (Mínimos Cuadrados)')
        ax2.set_xlabel("Tiempo (s)")
        ax2.set_ylabel("Temperatura (°C)")
        ax2.set_title("Curva de Calentamiento del Estator")
        ax2.grid(True, linestyle='--')
        ax2.legend()
        st.pyplot(fig2)

# ==========================================
# MÉTODO 3: NEWTON-RAPHSON
# ==========================================
with tab3:
    st.header("3. Predicción de Tiempo Crítico por Newton-Raphson")
    st.write("Búsqueda iterativa del instante $t$ en el cual la temperatura alcanza un umbral de control preventivo.")
    
    col3_1, col3_2 = st.columns([1, 2])
    
    with col3_1:
        st.subheader("Configuración de la Raíz")
        target_temp = st.slider("Umbral Térmico Preventivo (°C):", min_value=27.0, max_value=35.0, value=32.0, step=0.1)
        x0 = st.number_input("Valor Inicial de Semilla t0 (s):", value=50.0)
        
        t_curr = x0
        tol = 1e-6
        max_iter = 10
        tabla_nr = []
        
        for i in range(max_iter):
            f_val = (a0 + a1 * t_curr + a2 * (t_curr**2)) - target_temp
            f_prime = a1 + 2 * a2 * t_curr
            t_next = t_curr - (f_val / f_prime)
            
            tabla_nr.append({
                "Iteración": i + 1,
                "t_n (s)": round(t_curr, 6),
                "f(t_n)": round(f_val, 6),
                "f'(t_n)": round(f_prime, 6),
                "t_{n+1} (s)": round(t_next, 6),
                "Error": round(abs(t_next - t_curr), 6)
            })
            
            if abs(t_next - t_curr) < tol:
                break
            t_curr = t_next
            
        st.success(f"🎯 **Tiempo Crítico Encontrado:** t = **{t_next:.4f} s**")
        st.info("📌 Para 32.0 °C, el valor exacto del informe es **60.4622 s**.")

    with col3_2:
        st.subheader("Tabla de Convergencia Iterativa")
        df_nr = pd.DataFrame(tabla_nr)
        st.dataframe(df_nr, use_container_width=True)

# ==========================================
# MÉTODO 4: DERIVACIÓN NUMÉRICA
# ==========================================
with tab4:
    st.header("4. Tasa Instantánea de Calentamiento (Derivación Numérica)")
    st.write("Cálculo de la velocidad de variación de temperatura $dT/dt$ mediante Diferencias Finitas.")
    
    h_step = 5.0
    dT_dt = np.gradient(temp_exp, h_step)
    
    col4_1, col4_2 = st.columns([1, 2])
    
    with col4_1:
        st.subheader("Resultados Numéricos")
        df_der = pd.DataFrame({
            "Tiempo t (s)": t_exp,
            "Temp Te (°C)": temp_exp,
            "dT/dt (°C/s)": np.round(dT_dt, 4)
        })
        st.dataframe(df_der, use_container_width=True)
        
    with col4_2:
        st.subheader("Perfil de la Aceleración Térmica")
        fig4, ax4 = plt.subplots(figsize=(8, 4))
        ax4.plot(t_exp, dT_dt, 'm-o', linewidth=2, label='Tasa dT/dt (Diferencias Finitas)')
        ax4.set_xlabel("Tiempo t (s)")
        ax4.set_ylabel("Velocidad de Calentamiento (°C/s)")
        ax4.set_title("Respuesta Térmica Transitoria")
        ax4.grid(True, linestyle='--')
        ax4.legend()
        st.pyplot(fig4)

# ==========================================
# MÉTODO 5: INTEGRACIÓN Y EDO (RK4)
# ==========================================
with tab5:
    st.header("5. Simulación por Runge-Kutta 4to Orden (RK4) e Integración")
    
    col5_1, col5_2 = st.columns(2)
    
    with col5_1:
        st.subheader("A. Integración Numérica (Energía Total)")
        potencia = v_exp * i_exp  # P = V * I
        
        # Corrección de compatibilidad para versiones de NumPy
        try:
            energia_joules = float(np.trapezoid(potencia, t_exp))
        except AttributeError:
            energia_joules = float(np.trapz(potencia, t_exp))
            
        energia_wh = energia_joules / 3600.0
        
        st.metric("Energía Consumida (Joules)", f"{energia_joules:,.2f} J")
        st.metric("Energía Consumida (Watt-hora)", f"{energia_wh:.4f} Wh")
        st.caption("Calculado mediante la Regla del Trapecio sobre el perfil continuo de potencia.")

    with col5_2:
        st.subheader("B. Solución de EDO por RK4")
        tau = st.number_input("Constante de Tiempo Tau (s):", value=895.3)
        h_rk = 2.0
        
        t_rk = np.arange(0, 60 + h_rk, h_rk)
        v_rk = np.zeros(len(t_rk))
        v_rk[0] = 48.00
        
        for k in range(len(t_rk) - 1):
            f_rk = lambda v_val: -v_val / tau
            k1 = f_rk(v_rk[k])
            k2 = f_rk(v_rk[k] + 0.5 * h_rk * k1)
            k3 = f_rk(v_rk[k] + 0.5 * h_rk * k2)
            k4 = f_rk(v_rk[k] + h_rk * k3)
            v_rk[k+1] = v_rk[k] + (h_rk / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
            
        v_final_rk = v_rk[-1]
        st.metric("Voltaje Estimado a 60s (RK4)", f"{v_final_rk:.4f} V")
        st.write(f"Voltaje Real Medido a 60s: **44.7500 V**")
        st.write(f"Error Absoluto Marginal: **{abs(v_final_rk - 44.75):.4f} V**")

    st.subheader("Validación de la Dinámica de Voltaje (RK4 vs. Datos Reales)")
    fig5, ax5 = plt.subplots(figsize=(10, 4))
    ax5.plot(t_exp, v_exp, 'ro', markersize=8, label='Datos Experimentales (Tabla 1)')
    ax5.plot(t_rk, v_rk, 'b-', linewidth=2, label='Modelo Simulado (RK4)')
    ax5.set_xlabel("Tiempo t (s)")
    ax5.set_ylabel("Voltaje Vdc (V)")
    ax5.grid(True, linestyle='--')
    ax5.legend()
    st.pyplot(fig5)
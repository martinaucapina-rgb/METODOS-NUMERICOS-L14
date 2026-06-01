# =================================================================
# DATOS DEL ENSAYO DE COMPRESIÓN
# =================================================================
V = [1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0]
P = [300.2, 242.1, 201.5, 169.8, 151.0, 134.3, 120.8, 107.9, 99.1, 93.4, 86.2, 79.5, 75.1]

h = 0.25
n = 12  # 12 subintervalos

print("=====================================================================")
print(" RESOLUCIÓN Y DESARROLLO MATEMÁTICO DE LOS MÉTODOS")
print("=====================================================================\n")

# 1. Regla del Trapecio Compuesta
suma_trapecio = P[0] + P[n]
for i in range(1, n):
    suma_trapecio += 2 * P[i]
W_trapecio = (h / 2.0) * suma_trapecio

print("1) REGLA DEL TRAPECIO COMPUESTA:")
print(f"Fórmula: W = (h/2) * [ P_0 + 2*(P_1 + P_2 + ... + P_11) + P_12 ]")
print(f"Sustitución: W = ({h}/2) * [ {P[0]} + 2*({sum(P[1:11])}) + {P[12]} ]")
print(f"Resultado: W = {W_trapecio:.4f} kJ\n")


# 2. Regla de Simpson 1/3 Compuesta
suma_simpson13 = P[0] + P[n]
for i in range(1, n):
    if i % 2 != 0:
        suma_simpson13 += 4 * P[i]
    else:
        suma_simpson13 += 2 * P[i]
W_simpson13 = (h / 3.0) * suma_simpson13

print("2) REGLA DE SIMPSON 1/3 COMPUESTA:")
print(f"Fórmula: W = (h/3) * [ P_0 + 4*(P_impares) + 2*(P_pares) + P_12 ]")
# Calculamos sumas rápidas solo para la impresión ordenada
sum_imp = sum(P[i] for i in range(1, n) if i % 2 != 0)
sum_par = sum(P[i] for i in range(1, n) if i % 2 == 0)
print(f"Sustitución: W = ({h}/3) * [ {P[0]} + 4*({sum_imp:.1f}) + 2*({sum_par:.1f}) + {P[12]} ]")
print(f"Resultado: W = {W_simpson13:.4f} kJ\n")


# 3. Regla de Simpson 3/8 Compuesta
suma_simpson38 = P[0] + P[n]
for i in range(1, n):
    if i % 3 == 0:
        suma_simpson38 += 2 * P[i]
    else:
        suma_simpson38 += 3 * P[i]
W_simpson38 = ((3.0 * h) / 8.0) * suma_simpson38

print("3) REGLA DE SIMPSON 3/8 COMPUESTA:")
print(f"Fórmula: W = (3h/8) * [ P_0 + 3*(P_no_múltiplos_3) + 2*(P_múltiplos_3) + P_12 ]")
sum_m3 = sum(P[i] for i in range(1, n) if i % 3 == 0)
sum_nm3 = sum(P[i] for i in range(1, n) if i % 3 != 0)
print(f"Sustitución: W = (3*{h}/8) * [ {P[0]} + 3*({sum_nm3:.1f}) + 2*({sum_m3:.1f}) + {P[12]} ]")
print(f"Resultado: W = {W_simpson38:.4f} kJ\n")


# =================================================================
# INTERPOLACIÓN POR SPLINE CÚBICO NATURAL (PARA GAUSS-LEGENDRE)
# =================================================================
def thomas_tridiagonal(A_mat, B_mat, C_mat, D_mat):
    size = len(D_mat)
    cp, dp = [0.0] * size, [0.0] * size
    cp[0] = C_mat[0] / B_mat[0]
    dp[0] = D_mat[0] / B_mat[0]
    for i in range(1, size):
        den = B_mat[i] - A_mat[i] * cp[i-1]
        if i < size - 1: 
            cp[i] = C_mat[i] / den
        dp[i] = (D_mat[i] - A_mat[i] * dp[i-1]) / den
    x_sol = [0.0] * size
    x_sol[-1] = dp[-1]
    for i in range(size-2, -1, -1):
        x_sol[i] = dp[i] - cp[i] * x_sol[i+1]
    return x_sol

A_m, B_m, C_m, D_m = [0.0]*(n+1), [1.0]*(n+1), [0.0]*(n+1), [0.0]*(n+1)
for i in range(1, n):
    A_m[i] = h / 6.0
    B_m[i] = (2.0 * h) / 3.0
    C_m[i] = h / 6.0
    D_m[i] = ((P[i+1] - P[i]) / h) - ((P[i] - P[i-1]) / h)

M = thomas_tridiagonal(A_m, B_m, C_m, D_m)

def P_continua(V_obj):
    for i in range(n):
        if V[i] <= V_obj <= V[i+1]:
            t1 = M[i] * ((V[i+1] - V_obj)**3) / (6.0 * h)
            t2 = M[i+1] * ((V_obj - V[i])**3) / (6.0 * h)
            t3 = (P[i] / h - M[i] * h / 6.0) * (V[i+1] - V_obj)
            t4 = (P[i+1] / h - M[i+1] * h / 6.0) * (V_obj - V[i])
            return t1 + t2 + t3 + t4
    return 0.0

# =================================================================
# CUADRATURA DE GAUSS-LEGENDRE (m = 3 PUNTOS)
# =================================================================
a_g, b_g = 1.0, 4.0 
t_gauss = [-0.7745966692414834, 0.0, 0.7745966692414834]
w_gauss = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]

print("4) CUADRATURA DE GAUSS-LEGENDRE (m = 3 PUNTOS):")
print(f"Cambio de Variable V_i = ((4.0 - 1.0)/2)*t_i + ((4.0 + 1.0)/2)")
W_gauss = 0.0
for i in range(3):
    V_i = ((b_g - a_g) / 2.0) * t_gauss[i] + ((b_g + a_g) / 2.0)
    P_i = P_continua(V_i)
    W_gauss += w_gauss[i] * P_i
    print(f"  > Eval {i+1}: t = {t_gauss[i]:.4f} -> V = {V_i:.4f} L -> P(V) = {P_i:.2f} kPa")
W_gauss = ((b_g - a_g) / 2.0) * W_gauss
print(f"Resultado: W = (1.5) * [{w_gauss[0]:.4f}*P(V_1) + {w_gauss[1]:.4f}*P(V_2) + {w_gauss[2]:.4f}*P(V_3)] = {W_gauss:.4f} kJ\n")


# =================================================================
# CÁLCULO DE ERRORES RELATIVOS % (Tomando Gauss como referencia)
# =================================================================
err_trap = abs(W_gauss - W_trapecio) / W_gauss * 100
err_s13  = abs(W_gauss - W_simpson13) / W_gauss * 100
err_s38  = abs(W_gauss - W_simpson38) / W_gauss * 100
err_gauss = 0.0

# =================================================================
# IMPRESIÓN DE LA TABLA COMPARATIVA FORMATEADA
# =================================================================
print("="*60)
print("TABLA COMPARATIVA DE RESULTADOS")
print("="*60)
print("MÉTODO           | VALOR APROXIMADO | ERROR RELATIVOS %")
print("------------------------------------------------------------")
print(f"Trapecio         |     {W_trapecio:.4f}       |     {err_trap:.4f}%")
print(f"Simpson 1/3      |     {W_simpson13:.4f}       |     {err_s13:.4f}%")
print(f"Simpson 3/8      |     {W_simpson38:.4f}       |     {err_s38:.4f}%")
print(f"Gauss-Legendre   |     {W_gauss:.4f}       |     {err_gauss:.4f}%")
print("------------------------------------------------------------")
f_data = [
    10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5,
    35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5,
    60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5,
    85.0, 87.5, 90.0, 92.5, 95.0, 97.5, 100.0, 102.5, 105.0, 107.5
]

V_data = [
    0.842, 0.911, 0.986, 1.062, 1.143, 1.227, 1.314, 1.401, 1.482, 1.551,
    1.216, 1.048, 0.866, 0.689, 0.521, 0.364, 0.223, 0.103, 0.012, -0.041,
    -0.057, -0.034, 0.018, 0.096, 0.197, 0.318, 0.452, 0.579, 0.700, 0.809,
    0.611, 0.688, 0.756, 0.811, 0.856, 0.894, 0.926, 0.954, 0.980, 1.004
]

h = 2.5  

# =================================================================
# 1. VALORES PREVIOS DE DIFERENCIAS FINITAS (IND. 1 Y 2)
# =================================================================
# Extremo 10.0 kHz (Progresiva O(h²))
prog2_10 = (-3 * V_data[0] + 4 * V_data[1] - V_data[2]) / (2 * h)

# Punto 40.0 kHz
idx_40 = f_data.index(40.0)
dc2_40 = (V_data[idx_40 + 1] - V_data[idx_40 - 1]) / (2 * h)
dc4_40 = (-V_data[idx_40 + 2] + 8 * V_data[idx_40 + 1] - 8 * V_data[idx_40 - 1] + V_data[idx_40 - 2]) / (12 * h)

# Punto 70.0 kHz
idx_70 = f_data.index(70.0)
dc2_70 = (V_data[idx_70 + 1] - V_data[idx_70 - 1]) / (2 * h)
dc4_70 = (-V_data[idx_70 + 2] + 8 * V_data[idx_70 + 1] - 8 * V_data[idx_70 - 1] + V_data[idx_70 - 2]) / (12 * h)

# Punto 100.0 kHz
idx_100 = f_data.index(100.0)
dc2_100 = (V_data[idx_100 + 1] - V_data[idx_100 - 1]) / (2 * h)
dc4_100 = (-V_data[idx_100 + 2] + 8 * V_data[idx_100 + 1] - 8 * V_data[idx_100 - 1] + V_data[idx_100 - 2]) / (12 * h)


# =================================================================
# 4. MODELADO Y DERIVACIÓN DEL SPLINE CÚBICO NATURAL
# =================================================================
print("4. COMPARATIVA DE DIFERENCIAS FINITAS VS SPLINE CÚBICO")
print("=====================================================================")

def thomas_tridiagonal(A, B, C, D):
    n = len(D)
    cp, dp = [0.0] * n, [0.0] * n
    cp[0] = C[0] / B[0]
    dp[0] = D[0] / B[0]
    for i in range(1, n):
        den = B[i] - A[i] * cp[i-1]
        if i < n - 1: cp[i] = C[i] / den
        dp[i] = (D[i] - A[i] * dp[i-1]) / den
    x = [0.0] * n
    x[-1] = dp[-1]
    for i in range(n-2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i+1]
    return x

n_spl = len(f_data) - 1
A_m, B_m, C_m, D_m = [0.0]*(n_spl+1), [1.0]*(n_spl+1), [0.0]*(n_spl+1), [0.0]*(n_spl+1)

for i in range(1, n_spl):
    A_m[i] = h / 6.0
    B_m[i] = (2 * h) / 3.0
    C_m[i] = h / 6.0
    D_m[i] = ((V_data[i+1] - V_data[i]) / h) - ((V_data[i] - V_data[i-1]) / h)
    
M_V = thomas_tridiagonal(A_m, B_m, C_m, D_m)

def derivada_spline(x_obj, x, y, M):
    for i in range(len(x) - 1):
        if x[i] <= x_obj <= x[i+1]:
            t1 = -M[i] * ((x[i+1] - x_obj)**2) / (2 * h)
            t2 = M[i+1] * ((x_obj - x[i])**2) / (2 * h)
            t3 = (y[i+1] - y[i]) / h - (M[i+1] - M[i]) * h / 6.0
            return t1 + t2 + t3
    return 0.0

spl_10 = derivada_spline(10.0, f_data, V_data, M_V)
spl_40 = derivada_spline(40.0, f_data, V_data, M_V)
spl_70 = derivada_spline(70.0, f_data, V_data, M_V)
spl_100 = derivada_spline(100.0, f_data, V_data, M_V)

print("TABLA COMPARATIVA FINAL DE DERIVADAS (dV/df)")
print("+--------------+------------------+------------------+------------------+")
print("| Frecuencia   | Diferencia O(h²) |  Centrada O(h⁴)  |  Derivada Spline |")
print("+--------------+------------------+------------------+------------------+")
print("|   10.0 kHz   |     {:>8.4f}     |  N/A (Extremo)   |     {:>8.4f}     |".format(prog2_10, spl_10))
print("|   40.0 kHz   |     {:>8.4f}     |     {:>8.4f}     |     {:>8.4f}     |".format(dc2_40, dc4_40, spl_40))
print("|   70.0 kHz   |      {:>7.4f}     |      {:>7.4f}     |      {:>7.4f}     |".format(dc2_70, dc4_70, spl_70))
print("|  100.0 kHz   |      {:>7.4f}     |      {:>7.4f}     |      {:>7.4f}     |".format(dc2_100, dc4_100, spl_100))
print("+--------------+------------------+------------------+------------------+")
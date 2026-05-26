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

Z_data = [
    182.4, 178.9, 175.1, 171.0, 166.8, 162.7, 158.9, 155.4, 152.0, 149.0,
    146.1, 145.2, 145.8, 147.3, 149.9, 153.5, 158.0, 163.2, 168.9, 174.8,
    180.5, 186.2, 191.5, 196.2, 200.1, 203.1, 205.2, 206.3, 206.1, 204.7,
    198.0, 194.4, 190.9, 187.8, 185.1, 183.0, 181.6, 180.8, 180.6, 180.9
]

def resolver_lagrange(f_obj, x_data, y_data, var_letter):
    distancias = [abs(x - f_obj) for x in x_data]
    indices = sorted(range(len(distancias)), key=lambda k: distancias[k])[:3]
    indices.sort()
    
    f0, f1, f2 = x_data[indices[0]], x_data[indices[1]], x_data[indices[2]]
    y0, y1, y2 = y_data[indices[0]], y_data[indices[1]], y_data[indices[2]]
    
    L0 = ((f_obj - f1) * (f_obj - f2)) / ((f0 - f1) * (f0 - f2))
    L1 = ((f_obj - f0) * (f_obj - f2)) / ((f1 - f0) * (f1 - f2))
    L2 = ((f_obj - f0) * (f_obj - f1)) / ((f2 - f0) * (f2 - f1))
    
    valor_estimado = (L0 * y0) + (L1 * y1) + (L2 * y2)
    
    print(f"  > Nodos [{f0}, {f1}, {f2}] -> L0={L0:.4f}, L1={L1:.4f}, L2={L2:.4f} -> {var_letter}({f_obj} kHz) = {valor_estimado:.4f}")
    return valor_estimado

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

def resolver_spline(x, y, x_obj, var_letter):
    n = len(x) - 1
    h = [x[i+1] - x[i] for i in range(n)]
    A, B, C, D = [0.0]*(n+1), [1.0]*(n+1), [0.0]*(n+1), [0.0]*(n+1)
    
    for i in range(1, n):
        A[i] = h[i-1] / 6.0
        B[i] = (h[i-1] + h[i]) / 3.0
        C[i] = h[i] / 6.0
        D[i] = ((y[i+1] - y[i]) / h[i]) - ((y[i] - y[i-1]) / h[i-1])
        
    M = thomas_tridiagonal(A, B, C, D)
    
    tramo = 0
    for i in range(n):
        if x[i] <= x_obj <= x[i+1]:
            tramo = i
            break
            
    i = tramo
    hi = h[i]
    t1 = M[i] * ((x[i+1] - x_obj)**3) / (6.0 * hi)
    t2 = M[i+1] * ((x_obj - x[i])**3) / (6.0 * hi)
    t3 = (y[i] - (M[i] * (hi**2)) / 6.0) * ((x[i+1] - x_obj) / hi)
    t4 = (y[i+1] - (M[i+1] * (hi**2)) / 6.0) * ((x_obj - x[i]) / hi)
    valor_estimado = t1 + t2 + t3 + t4
    
    print(f"  > Tramo {i} [{x[i]} a {x[i+1]}] -> M[{i}]={M[i]:.4f}, M[{i+1}]={M[i+1]:.4f} -> {var_letter}({x_obj} kHz) = {valor_estimado:.4f}")
    return valor_estimado

# 1. BLOQUE DE INTERPOLACIÓN DE LAGRANGE
print("=================================================================")
print("PROCESO 1: INTERPOLACIÓN DE LAGRANGE DE SEGUNDO GRADO (LOCAL)")
print("=================================================================")
v_lag_41 = resolver_lagrange(41.0, f_data, V_data, "V")
z_lag_41 = resolver_lagrange(41.0, f_data, Z_data, "Z")
v_lag_73 = resolver_lagrange(73.0, f_data, V_data, "V")
z_lag_73 = resolver_lagrange(73.0, f_data, Z_data, "Z")
print()

# 2. BLOQUE DE SPLINE CÚBICO NATURAL
print("=================================================================")
print("PROCESO 2: MODELADO GLOBAL CON SPLINE CÚBICO NATURAL")
print("=================================================================")
v_spl_41 = resolver_spline(f_data, V_data, 41.0, "V")
z_spl_41 = resolver_spline(f_data, Z_data, 41.0, "Z")
v_spl_73 = resolver_spline(f_data, V_data, 73.0, "V")
z_spl_73 = resolver_spline(f_data, Z_data, 73.0, "Z")
print()

# 3. TABLA COMPARATIVA DE RESULTADOS FINALES
print("TABLA COMPARATIVA DE RESULTADOS DE INTERPOLACIÓN")
print("+" + "-"*71 + "+")
print("| {:^18} | {:^22} | {:^23} |".format("Variable (kHz)", "Lagrange (Grado 2)", "Spline Cúbico Natural"))
print("+" + "-"*71 + "+")
print("| {:<18} | {:^22.4f} | {:^23.4f} |".format("V(41.0 kHz)", v_lag_41, v_spl_41))
print("| {:<18} | {:^22.4f} | {:^23.4f} |".format("Z(41.0 kHz)", z_lag_41, z_spl_41))
print("+" + "-"*71 + "+")
print("| {:<18} | {:^22.4f} | {:^23.4f} |".format("V(73.0 kHz)", v_lag_73, v_spl_73))
print("| {:<18} | {:^22.4f} | {:^23.4f} |".format("Z(73.0 kHz)", z_lag_73, z_spl_73))
print("+" + "-"*71 + "+")
print()